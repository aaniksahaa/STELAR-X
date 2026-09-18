/**
 * stelarx_similarity.cu
 * =====================
 * GPU similarity-matrix computation reproducing ASTRAL-MP's
 * SimilarityMatrix.populateByQuartetDistance byte-for-byte.
 *
 * FORMULA (bridge identity, validated 1141/1141 pairs):
 *   same_side_T(a, b)  =  C2(kt − 2)  −  QD_T(a, b)
 *
 * O(1) closed form for QD via Euler tour + RMQ:
 *   QD_T(x, y) = ½ · [ (F(x) − F(cx)) + (F(y) − F(cy))
 *                     + (cxS − 1)·Z + (cyS − 1)·Z ]
 * where w = LCA(x, y), cx = child of w on x-side, cy = child of w on y-side,
 *   cxS = s(cx),  cyS = s(cy),   Z = kt − cxS − cyS,
 *   F(v) the root→v path prefix described in EulerTourBuilder.
 *
 * O(1) LCA + child-payload query via Euler tour RMQ:
 *   fa = firstOcc[a],  fb = firstOcc[b];   l = min(fa,fb), r = max(...)
 *   k_lvl = floor(log2(r − l + 1));  l2 = r − 2^k_lvl + 1
 *   pL = sparseArgmin[k_lvl][l];  pR = sparseArgmin[k_lvl][l2]
 *   leftWins = (eulerDepth[pL] <= eulerDepth[pR])
 *   The selected position is the INTERMEDIATE visit of LCA(x,y). Child
 *   payloads s(LCA.left), F(LCA.left), s(LCA.right), F(LCA.right) are read
 *   from the base Euler arrays at that exact position. The leaf with the
 *   smaller firstOcc is in LCA.left.
 *
 * ARCHITECTURE:
 *   Δ-tree batching : tree data on GPU  O(Δ · n · log n)
 *   B×B pair tiling : output tile       O(B²) doubles
 *   No atomics: each thread owns a unique (a,b) cell.
 */

#include <cuda_runtime.h>
#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "stelarx_platform.h"

// ── Exact host output layout ─────────────────────────────────────────────────

struct SimHostOutput {
    JNIEnv* env = nullptr;
    bool packed = false;
    int segment_shift = 0;
    long long segment_mask = 0;
    jdoubleArray dense_num_ref = nullptr;
    jdoubleArray dense_den_ref = nullptr;
    jdouble* dense_num = nullptr;
    jdouble* dense_den = nullptr;
    std::vector<jdoubleArray> num_refs;
    std::vector<jdoubleArray> den_refs;
    std::vector<jdouble*> num_segments;
    std::vector<jdouble*> den_segments;

    ~SimHostOutput() { release(); }

    bool acquire(JNIEnv* e, jdoubleArray j_dense_num, jdoubleArray j_dense_den,
                 jobjectArray j_packed_num, jobjectArray j_packed_den, int shift) {
        env = e;
        packed = j_packed_num != nullptr;
        if (!packed) {
            dense_num_ref = j_dense_num;
            dense_den_ref = j_dense_den;
            jboolean is_copy;
            dense_num = env->GetDoubleArrayElements(j_dense_num, &is_copy);
            dense_den = env->GetDoubleArrayElements(j_dense_den, &is_copy);
            return dense_num != nullptr && dense_den != nullptr;
        }

        if (j_packed_den == nullptr || shift <= 0 || shift >= 62) {
            jclass ex = env->FindClass("java/lang/IllegalArgumentException");
            if (ex != nullptr) env->ThrowNew(ex, "invalid packed similarity output layout");
            return false;
        }
        segment_shift = shift;
        segment_mask = (1LL << shift) - 1LL;
        jsize count = env->GetArrayLength(j_packed_num);
        if (env->GetArrayLength(j_packed_den) != count) {
            jclass ex = env->FindClass("java/lang/IllegalArgumentException");
            if (ex != nullptr) env->ThrowNew(ex, "packed similarity segment counts differ");
            return false;
        }
        num_refs.reserve(count); den_refs.reserve(count);
        num_segments.reserve(count); den_segments.reserve(count);
        for (jsize s = 0; s < count; s++) {
            jdoubleArray nr = (jdoubleArray)env->GetObjectArrayElement(j_packed_num, s);
            jdoubleArray dr = (jdoubleArray)env->GetObjectArrayElement(j_packed_den, s);
            jboolean is_copy;
            jdouble* np = env->GetDoubleArrayElements(nr, &is_copy);
            jdouble* dp = env->GetDoubleArrayElements(dr, &is_copy);
            if (np == nullptr || dp == nullptr) {
                if (np != nullptr) env->ReleaseDoubleArrayElements(nr, np, 0);
                if (dp != nullptr) env->ReleaseDoubleArrayElements(dr, dp, 0);
                if (nr != nullptr) env->DeleteLocalRef(nr);
                if (dr != nullptr) env->DeleteLocalRef(dr);
                return false;
            }
            num_refs.push_back(nr); den_refs.push_back(dr);
            num_segments.push_back(np); den_segments.push_back(dp);
        }
        return true;
    }

    __host__ inline void addPacked(int a, int b, int n, double nv, double dv) {
        int i = a, j = b;
        if (i > j) { int t = i; i = j; j = t; }
        long long index = (long long)i * n - (long long)i * (i + 1LL) / 2LL + j;
        int segment = (int)(index >> segment_shift);
        long long offset = index & segment_mask;
        num_segments[segment][offset] += nv;
        den_segments[segment][offset] += dv;
    }

    void release() {
        if (!env) return;
        if (!packed) {
            if (dense_num) env->ReleaseDoubleArrayElements(dense_num_ref, dense_num, 0);
            if (dense_den) env->ReleaseDoubleArrayElements(dense_den_ref, dense_den, 0);
            dense_num = dense_den = nullptr;
            return;
        }
        for (size_t s = 0; s < num_segments.size(); s++) {
            env->ReleaseDoubleArrayElements(num_refs[s], num_segments[s], 0);
            env->ReleaseDoubleArrayElements(den_refs[s], den_segments[s], 0);
            env->DeleteLocalRef(num_refs[s]);
            env->DeleteLocalRef(den_refs[s]);
        }
        num_segments.clear(); den_segments.clear();
        num_refs.clear(); den_refs.clear();
    }
};

// ── Utility: timing ──────────────────────────────────────────────────────────

static double sim_now_sec() {
    return stelarx_now_sec();
}

static void sim_fmt_duration(double sec, char* buf, int bufsz) {
    int s = (int)sec, m = s / 60;
    s %= 60;
    if (m >= 60) snprintf(buf, bufsz, "%dh%02dm%02ds", m/60, m%60, s);
    else         snprintf(buf, bufsz, "%02d:%02d",     m, s);
}

static int sim_use_color() {
    if (getenv("NO_COLOR"))    return 0;
    if (getenv("FORCE_COLOR")) return 1;
    return stelarx_stderr_isatty();
}

static void sim_fmt_rate(double rate, char* buf, int bufsz) {
    if (rate <= 0)       snprintf(buf, bufsz, "?it/s");
    else if (rate >= 1)  snprintf(buf, bufsz, "%.1fit/s", rate);
    else                 snprintf(buf, bufsz, "%.2fs/it", 1.0 / rate);
}

static void sim_print_progress(int work_done, int total_work, double elapsed,
                                int color, int is_last) {
    double pct  = (total_work > 0) ? 100.0 * work_done / total_work : 100.0;
    double rate = (work_done > 0 && elapsed > 0) ? work_done / elapsed : 0.0;
    double eta  = (rate > 0 && work_done < total_work) ? (total_work - work_done) / rate : 0.0;

    char elapsed_buf[32], eta_buf[32], rate_buf[32];
    sim_fmt_duration(elapsed, elapsed_buf, sizeof(elapsed_buf));
    sim_fmt_duration(eta,     eta_buf,     sizeof(eta_buf));
    sim_fmt_rate(rate, rate_buf, sizeof(rate_buf));

    const int BAR_W = 28;
    int filled = (int)(BAR_W * pct / 100.0);
    const char* FULL  = "\xE2\x96\x88";
    const char* EMPTY = "\xE2\x96\x91";
    char bar[4 + BAR_W * 3 + 4];
    int pos = 0;
    bar[pos++] = '[';
    for (int i = 0; i < BAR_W; i++) {
        const char* ch = (i < filled) ? FULL : EMPTY;
        bar[pos++] = ch[0]; bar[pos++] = ch[1]; bar[pos++] = ch[2];
    }
    bar[pos++] = ']'; bar[pos] = '\0';

    if (color)
        fprintf(stderr,
            "     \033[2m▸  \033[0mSimilarity matrix (GPU)  "
            "\033[32m%s\033[0m  %d/%d (%d%%)"
            "  \033[2m[%s<%s, \033[0m\033[33m%s\033[0m\033[2m]\033[0m\r",
            bar, work_done, total_work, (int)pct,
            elapsed_buf, work_done > 0 ? eta_buf : "?", rate_buf);
    else
        fprintf(stderr,
            "     ▸  Similarity matrix (GPU)  %s  %d/%d (%d%%)  [%s<%s, %s]\r",
            bar, work_done, total_work, (int)pct,
            elapsed_buf, work_done > 0 ? eta_buf : "?", rate_buf);
    fflush(stderr);
    if (is_last) fprintf(stderr, "\n");
}

// ── Utility: CUDA error check ─────────────────────────────────────────────────

#define SIM_CUDA_CHECK(call) \
    do { \
        cudaError_t _e = (call); \
        if (_e != cudaSuccess) { \
            char _msg[768]; \
            snprintf(_msg, sizeof(_msg), \
                     "STELAR-X similarity CUDA error in %s at %s:%d: %s", \
                     #call, __FILE__, __LINE__, cudaGetErrorString(_e)); \
            fprintf(stderr, "[STELAR-X sim] %s\n", _msg); \
            jclass _ex = env->FindClass("java/lang/RuntimeException"); \
            if (_ex != nullptr) env->ThrowNew(_ex, _msg); \
            return; \
        } \
    } while(0)

// ── GPU Kernel ────────────────────────────────────────────────────────────────

/**
 * sim_tile_kernel
 * ───────────────
 * One thread per (da, db) = position within the B×B pair tile.
 * Global taxa: a = a0+da, b = b0+db.
 *
 * For each tree t in the current batch:
 *   1. Presence check via firstOcc.
 *   2. Skip if kt < 4 (C2(kt-2) = 0).
 *   3. O(1) RMQ → LCA's intermediate-position child payloads.
 *   4. Compute twoQD = 2·QD via the closed form.
 *   5. ss = C2(kt-2) − twoQD/2;  numAcc += ss;  denAcc += C2(kt-2).
 */
__global__ void sim_tile_kernel(
    const short*  __restrict__ euler_depths,        // [delta * E_max]
    const double* __restrict__ euler_F,             // [delta * E_max]
    const short*  __restrict__ euler_left_child_s,  // [delta * E_max]
    const double* __restrict__ euler_left_child_f,  // [delta * E_max]
    const short*  __restrict__ euler_right_child_s, // [delta * E_max]
    const double* __restrict__ euler_right_child_f, // [delta * E_max]
    const unsigned short* __restrict__ sparse_argmin, // [delta * LOG * E_max]
    const int*    __restrict__ first_occ,           // [delta * n]
    const int*    __restrict__ leaf_count,          // [delta]
    int delta, int n, int E_max, int LOG,
    int a0, int b0, int bA, int bB,
    double* __restrict__ tile_num,
    double* __restrict__ tile_den
) {
    int da = blockIdx.x * blockDim.x + threadIdx.x;
    int db = blockIdx.y * blockDim.y + threadIdx.y;
    if (da >= bA || db >= bB) return;

    int a = a0 + da;
    int b = b0 + db;
    if (a >= n || b >= n || a >= b) return;   // upper triangle only

    double local_num = 0.0;
    double local_den = 0.0;

    for (int t = 0; t < delta; t++) {
        long focc_off = (long)t * n;
        int fa = first_occ[focc_off + a];
        int fb = first_occ[focc_off + b];
        if (fa < 0 || fb < 0) continue;

        int kt = leaf_count[t];
        long long cc = (long long)(kt - 2) * (kt - 3) / 2;   // C2(kt-2)
        if (cc <= 0) continue;

        // ── O(1) RMQ: select left-biased argmin position ────────────────────
        int l   = (fa < fb) ? fa : fb;
        int r   = (fa < fb) ? fb : fa;
        int len = r - l + 1;
        int k_lvl = 31 - __clz(len);
        int l2    = r - (1 << k_lvl) + 1;

        long sp_off = (long)t * LOG * E_max;
        long ol     = sp_off + (long)k_lvl * E_max + l;
        long ol2    = sp_off + (long)k_lvl * E_max + l2;

        // Each sparse cell stores only its interval's left-biased argmin Euler
        // position. Compare the two overlap candidates exactly as before, then
        // fetch the winning child payloads from the base Euler arrays.
        int posL = (int)sparse_argmin[ol];
        int posR = (int)sparse_argmin[ol2];
        long ed_off = (long)t * E_max;
        short dL = euler_depths[ed_off + posL];
        short dR = euler_depths[ed_off + posR];
        int pickPos = (dL <= dR) ? posL : posR;  // identical left-biased tie rule
        long pickIdx = ed_off + pickPos;

        int    leftS = (int)   euler_left_child_s [pickIdx];
        double leftF = (double)euler_left_child_f [pickIdx];
        int    rightS= (int)   euler_right_child_s[pickIdx];
        double rightF= (double)euler_right_child_f[pickIdx];

        // Map (leftLeaf, rightLeaf) by tour order back to (a, b).
        // The leaf with the smaller firstOcc is in the LCA's LEFT child.
        int    aS;   double aF;
        int    bS;   double bF;
        if (fa <= fb) {
            aS = leftS;   aF = leftF;
            bS = rightS;  bF = rightF;
        } else {
            aS = rightS;  aF = rightF;
            bS = leftS;   bF = leftF;
        }

        double Fa = euler_F[ed_off + fa];
        double Fb = euler_F[ed_off + fb];

        long long Z = (long long)(kt - aS - bS);
        // twoQD = (Fa - aF) + (Fb - bF) + (aS - 1)*Z + (bS - 1)*Z
        double twoQD = (Fa - aF) + (Fb - bF)
                     + (double)((long long)(aS - 1) * Z)
                     + (double)((long long)(bS - 1) * Z);

        double ss = (double)cc - twoQD * 0.5;     // same-side count
        local_num += ss;
        local_den += (double)cc;
    }

    long tidx = (long)da * bB + db;
    tile_num[tidx] += local_num;
    tile_den[tidx] += local_den;
}

// ── Wide blocked-RMQ GPU kernel ──────────────────────────────────────────────

__device__ __forceinline__ int sim_wide_micro_argmin(
    const int* __restrict__ depths,
    const unsigned char* __restrict__ micro,
    long long ed_off, long long micro_off,
    int E_max, int block_size,
    int lo, int hi
) {
    int width = hi - lo + 1;
    int lvl = 31 - __clz(width);
    int second = hi - (1 << lvl) + 1;
    int block_start = (lo / block_size) * block_size;
    long long row = micro_off + (long long)lvl * E_max;
    int pos_l = block_start + (int)micro[row + lo];
    int pos_r = block_start + (int)micro[row + second];
    return (depths[ed_off + pos_l] <= depths[ed_off + pos_r]) ? pos_l : pos_r;
}

__device__ __forceinline__ int sim_wide_argmin(
    const int* __restrict__ depths,
    const unsigned char* __restrict__ micro,
    const int* __restrict__ macro,
    int tree, int E_max, int micro_log,
    int block_size, int block_max, int macro_log,
    int lo, int hi
) {
    long long ed_off = (long long)tree * E_max;
    long long micro_off = (long long)tree * micro_log * E_max;
    int left_block = lo / block_size;
    int right_block = hi / block_size;
    if (left_block == right_block) {
        return sim_wide_micro_argmin(depths, micro, ed_off, micro_off,
            E_max, block_size, lo, hi);
    }

    int best = sim_wide_micro_argmin(depths, micro, ed_off, micro_off,
        E_max, block_size, lo, (left_block + 1) * block_size - 1);

    int first_whole = left_block + 1;
    int last_whole = right_block - 1;
    if (first_whole <= last_whole) {
        int count = last_whole - first_whole + 1;
        int lvl = 31 - __clz(count);
        int second = last_whole - (1 << lvl) + 1;
        long long macro_off = (long long)tree * macro_log * block_max
                            + (long long)lvl * block_max;
        int pos_l = macro[macro_off + first_whole];
        int pos_r = macro[macro_off + second];
        int middle = (depths[ed_off + pos_l] <= depths[ed_off + pos_r])
                   ? pos_l : pos_r;
        if (depths[ed_off + middle] < depths[ed_off + best]) best = middle;
    }

    int right = sim_wide_micro_argmin(depths, micro, ed_off, micro_off,
        E_max, block_size, right_block * block_size, hi);
    if (depths[ed_off + right] < depths[ed_off + best]) best = right;
    return best;
}

__global__ void sim_tile_kernel_wide(
    const int*    __restrict__ euler_depths,
    const double* __restrict__ euler_F,
    const int*    __restrict__ euler_left_child_s,
    const double* __restrict__ euler_left_child_f,
    const int*    __restrict__ euler_right_child_s,
    const double* __restrict__ euler_right_child_f,
    const unsigned char* __restrict__ micro_argmin,
    const int*    __restrict__ macro_argmin,
    const int*    __restrict__ first_occ,
    const int*    __restrict__ leaf_count,
    int delta, int n, int E_max, int micro_log,
    int block_size, int block_max, int macro_log,
    int a0, int b0, int bA, int bB,
    double* __restrict__ tile_num,
    double* __restrict__ tile_den
) {
    int da = blockIdx.x * blockDim.x + threadIdx.x;
    int db = blockIdx.y * blockDim.y + threadIdx.y;
    if (da >= bA || db >= bB) return;

    int a = a0 + da;
    int b = b0 + db;
    if (a >= n || b >= n || a >= b) return;

    double local_num = 0.0;
    double local_den = 0.0;
    for (int t = 0; t < delta; t++) {
        long long focc_off = (long long)t * n;
        int fa = first_occ[focc_off + a];
        int fb = first_occ[focc_off + b];
        if (fa < 0 || fb < 0) continue;

        int kt = leaf_count[t];
        long long cc = (long long)(kt - 2) * (kt - 3) / 2;
        if (cc <= 0) continue;

        int lo = (fa < fb) ? fa : fb;
        int hi = (fa < fb) ? fb : fa;
        int pick_pos = sim_wide_argmin(euler_depths, micro_argmin, macro_argmin,
            t, E_max, micro_log, block_size, block_max, macro_log, lo, hi);
        long long ed_off = (long long)t * E_max;
        long long pick_idx = ed_off + pick_pos;

        int leftS = euler_left_child_s[pick_idx];
        double leftF = euler_left_child_f[pick_idx];
        int rightS = euler_right_child_s[pick_idx];
        double rightF = euler_right_child_f[pick_idx];

        int aS, bS;
        double aF, bF;
        if (fa <= fb) {
            aS = leftS; aF = leftF;
            bS = rightS; bF = rightF;
        } else {
            aS = rightS; aF = rightF;
            bS = leftS; bF = leftF;
        }

        double Fa = euler_F[ed_off + fa];
        double Fb = euler_F[ed_off + fb];
        long long Z = (long long)(kt - aS - bS);
        double twoQD = (Fa - aF) + (Fb - bF)
                     + (double)((long long)(aS - 1) * Z)
                     + (double)((long long)(bS - 1) * Z);
        local_num += (double)cc - twoQD * 0.5;
        local_den += (double)cc;
    }

    long tidx = (long)da * bB + db;
    tile_num[tidx] += local_num;
    tile_den[tidx] += local_den;
}

// ── JNI entry point ───────────────────────────────────────────────────────────

extern "C" JNIEXPORT void JNICALL
Java_stelarx_gpu_GPUSimilarityMatrix_computeSimilarityGPU(
    JNIEnv*  env,
    jclass   cls,
    jshortArray  j_euler_depths,
    jdoubleArray j_euler_F,
    jshortArray  j_euler_left_child_s,
    jdoubleArray j_euler_left_child_f,
    jshortArray  j_euler_right_child_s,
    jdoubleArray j_euler_right_child_f,
    jcharArray   j_sparse_argmin,
    jintArray    j_first_occ,
    jintArray    j_euler_len,
    jintArray    j_leaf_count,
    jint     numTrees,
    jint     n,
    jint     E_max,
    jint     LOG,
    jint     tileSizeB,
    jint     treeVramCapMiB,
    jdouble  progressInterval,
    jint     progressMaxSteps,
    jdoubleArray j_num_sum_out,
    jdoubleArray j_den_sum_out,
    jobjectArray j_packed_num_out,
    jobjectArray j_packed_den_out,
    jint packed_segment_shift
) {
    jboolean isCopy;
    jshort*  h_euler   = env->GetShortArrayElements (j_euler_depths,         &isCopy);
    jdouble* h_eulerF  = env->GetDoubleArrayElements(j_euler_F,              &isCopy);
    jshort*  h_eLcS    = env->GetShortArrayElements (j_euler_left_child_s,   &isCopy);
    jdouble* h_eLcF    = env->GetDoubleArrayElements(j_euler_left_child_f,   &isCopy);
    jshort*  h_eRcS    = env->GetShortArrayElements (j_euler_right_child_s,  &isCopy);
    jdouble* h_eRcF    = env->GetDoubleArrayElements(j_euler_right_child_f,  &isCopy);
    jchar*   h_argmin  = env->GetCharArrayElements  (j_sparse_argmin,        &isCopy);
    jint*    h_focc    = env->GetIntArrayElements   (j_first_occ,            &isCopy);
    jint*    h_elen    = env->GetIntArrayElements   (j_euler_len,            &isCopy);
    jint*    h_lcount  = env->GetIntArrayElements   (j_leaf_count,           &isCopy);
    SimHostOutput output;
    if (!output.acquire(env, j_num_sum_out, j_den_sum_out,
                        j_packed_num_out, j_packed_den_out, packed_segment_shift)) return;

    size_t free_vram = 0, total_vram = 0;
    cudaMemGetInfo(&free_vram, &total_vram);

    int B = tileSizeB;
    if (B <= 0) {
        B = (int)ceil(sqrt((double)n * numTrees));
        if (B > n) B = n;
    }
    while (B > 1 && 2LL * B * B * (long long)sizeof(double) > (long long)(free_vram * 0.40)) B /= 2;
    if (B < 1) B = 1;

    size_t tile_vram = 2ULL * B * B * sizeof(double);
    size_t headroom = 64ULL * 1024 * 1024;
    size_t available = (free_vram > tile_vram + headroom)
                     ? free_vram - tile_vram - headroom
                     : 16ULL * 1024 * 1024;
    size_t requested = (size_t)treeVramCapMiB * 1024 * 1024;
    size_t remaining = (requested < available) ? requested : available;

    // Per-tree bytes: see Java side comment in SimilarityMatrixBuilder.buildGPU.
    //   euler arrays:  E_max × (2 + 2 + 2 + 8 + 8 + 8) = E_max × 30
    //   sparse argmin: LOG × E_max × 2 (unsigned-16 Euler position)
    //   firstOcc:      n × 4
    //   leafCount:     4
    size_t per_tree = (size_t)E_max * 30
                    + (size_t)LOG * E_max * sizeof(unsigned short)
                    + (size_t)n * 4
                    + sizeof(int);
    int delta = (per_tree > 0) ? (int)(remaining / per_tree) : numTrees;
    if (delta < 1)        delta = 1;
    if (delta > numTrees) delta = numTrees;

    int num_batches    = (numTrees + delta - 1) / delta;
    int num_tiles_side = (n + B - 1) / B;
    int num_tiles      = num_tiles_side * (num_tiles_side + 1) / 2;

    fprintf(stderr,
        "\n[STELAR-X sim] GPU similarity matrix: n=%d  k=%d  "
        "tile B=%d  tree-batch Δ=%d  (%d batches × %d tiles)\n",
        n, numTrees, B, delta, num_batches, num_tiles);
    fprintf(stderr,
        "[STELAR-X sim] GPU VRAM: tile %.1f MB  tree-data %.1f MB  "
        "(cap %d MiB; free %.0f MB / total %.0f MB)\n",
        tile_vram / 1e6,
        (double)delta * per_tree / 1e6,
        treeVramCapMiB,
        free_vram / 1e6, total_vram / 1e6);
    if (num_batches > 1) {
        fprintf(stderr,
            "[STELAR-X sim] NOTE: similarity tree-data is bounded to %d MiB (%d batches). "
            "If this phase is a bottleneck, raise --gpu-sim-vram-cap-mb to use fewer batches; "
            "results are unchanged.\n",
            treeVramCapMiB, num_batches);
    }

    // ── Allocate GPU tree-data buffers ────────────────────────────────────────
    short*  d_euler   = nullptr;
    double* d_eulerF  = nullptr;
    short*  d_eLcS    = nullptr;
    double* d_eLcF    = nullptr;
    short*  d_eRcS    = nullptr;
    double* d_eRcF    = nullptr;
    unsigned short* d_argmin = nullptr;
    int*    d_focc    = nullptr;
    int*    d_lcount  = nullptr;

    long long sz_short_e  = (long long)delta * E_max * sizeof(short);
    long long sz_double_e = (long long)delta * E_max * sizeof(double);
    long long sz_argmin = (long long)delta * LOG * E_max * sizeof(unsigned short);

    SIM_CUDA_CHECK(cudaMalloc(&d_euler,  sz_short_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eulerF, sz_double_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eLcS,   sz_short_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eLcF,   sz_double_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eRcS,   sz_short_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eRcF,   sz_double_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_argmin, sz_argmin));
    SIM_CUDA_CHECK(cudaMalloc(&d_focc,   (long long)delta * n * sizeof(int)));
    SIM_CUDA_CHECK(cudaMalloc(&d_lcount, delta * sizeof(int)));

    double* d_tile_num = nullptr;
    double* d_tile_den = nullptr;
    SIM_CUDA_CHECK(cudaMalloc(&d_tile_num, (long long)B * B * sizeof(double)));
    SIM_CUDA_CHECK(cudaMalloc(&d_tile_den, (long long)B * B * sizeof(double)));

    double* h_tile_num = nullptr;
    double* h_tile_den = nullptr;
    SIM_CUDA_CHECK(cudaMallocHost(&h_tile_num, (long long)B * B * sizeof(double)));
    SIM_CUDA_CHECK(cudaMallocHost(&h_tile_den, (long long)B * B * sizeof(double)));

    int    total_work       = num_batches * num_tiles;
    int    work_done        = 0;
    double t_start          = sim_now_sec();
    double t_last_print     = t_start - progressInterval;
    double last_pct_printed = -1.0;
    const bool step_mode    = (progressMaxSteps > 0);
    int    use_color        = sim_use_color();

    for (int t0 = 0; t0 < numTrees; t0 += delta) {
        int dt = (t0 + delta > numTrees) ? numTrees - t0 : delta;

        long long off_e = (long long)t0 * E_max;
        long long bytes_short_e  = (long long)dt * E_max * sizeof(short);
        long long bytes_double_e = (long long)dt * E_max * sizeof(double);
        long long off_s = (long long)t0 * LOG * E_max;
        long long bytes_argmin = (long long)dt * LOG * E_max * sizeof(unsigned short);

        SIM_CUDA_CHECK(cudaMemcpy(d_euler,  h_euler  + off_e, bytes_short_e,  cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eulerF, h_eulerF + off_e, bytes_double_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eLcS,   h_eLcS   + off_e, bytes_short_e,  cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eLcF,   h_eLcF   + off_e, bytes_double_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eRcS,   h_eRcS   + off_e, bytes_short_e,  cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eRcF,   h_eRcF   + off_e, bytes_double_e, cudaMemcpyHostToDevice));

        SIM_CUDA_CHECK(cudaMemcpy(d_argmin, h_argmin + off_s, bytes_argmin, cudaMemcpyHostToDevice));

        SIM_CUDA_CHECK(cudaMemcpy(d_focc,
            h_focc + (long long)t0 * n,
            (long long)dt * n * sizeof(int), cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_lcount,
            h_lcount + t0,
            dt * sizeof(int), cudaMemcpyHostToDevice));

        for (int a0 = 0; a0 < n; a0 += B) {
            int bA = (a0 + B > n) ? n - a0 : B;
            for (int b0 = a0; b0 < n; b0 += B) {
                int bB_tile = (b0 + B > n) ? n - b0 : B;

                SIM_CUDA_CHECK(cudaMemset(d_tile_num, 0,
                    (long long)bA * bB_tile * sizeof(double)));
                SIM_CUDA_CHECK(cudaMemset(d_tile_den, 0,
                    (long long)bA * bB_tile * sizeof(double)));

                dim3 block(32, 32);
                dim3 grid((bA + 31) / 32, (bB_tile + 31) / 32);
                sim_tile_kernel<<<grid, block>>>(
                    d_euler, d_eulerF,
                    d_eLcS, d_eLcF, d_eRcS, d_eRcF,
                    d_argmin,
                    d_focc, d_lcount,
                    dt, n, E_max, LOG,
                    a0, b0, bA, bB_tile,
                    d_tile_num, d_tile_den
                );
                SIM_CUDA_CHECK(cudaDeviceSynchronize());

                SIM_CUDA_CHECK(cudaMemcpy(h_tile_num, d_tile_num,
                    (long long)bA * bB_tile * sizeof(double), cudaMemcpyDeviceToHost));
                SIM_CUDA_CHECK(cudaMemcpy(h_tile_den, d_tile_den,
                    (long long)bA * bB_tile * sizeof(double), cudaMemcpyDeviceToHost));

                if (!output.packed) {
                    for (int da = 0; da < bA; da++) {
                        int a = a0 + da;
                        for (int db = 0; db < bB_tile; db++) {
                            int b = b0 + db;
                            if (a >= b) continue;
                            long tidx = (long)da * bB_tile + db;
                            double nv = h_tile_num[tidx];
                            double dv = h_tile_den[tidx];
                            if (dv == 0.0) continue;
                            output.dense_num[a * n + b] += nv;
                            output.dense_num[b * n + a] += nv;
                            output.dense_den[a * n + b] += dv;
                            output.dense_den[b * n + a] += dv;
                        }
                    }
                } else {
                    for (int da = 0; da < bA; da++) {
                        int a = a0 + da;
                        for (int db = 0; db < bB_tile; db++) {
                            int b = b0 + db;
                            if (a >= b) continue;
                            long tidx = (long)da * bB_tile + db;
                            double nv = h_tile_num[tidx];
                            double dv = h_tile_den[tidx];
                            if (dv == 0.0) continue;
                            output.addPacked(a, b, n, nv, dv);
                        }
                    }
                }

                work_done++;
                double now = sim_now_sec();
                double pct = (total_work > 0) ? 100.0 * work_done / total_work : 100.0;
                bool is_last = (work_done == total_work);
                bool should_print = step_mode
                    ? (is_last || pct - last_pct_printed >= 100.0 / progressMaxSteps)
                    : (is_last || now - t_last_print >= progressInterval);
                if (should_print) {
                    sim_print_progress(work_done, total_work, now - t_start, use_color, is_last);
                    t_last_print     = now;
                    last_pct_printed = pct;
                }
            }
        }
    }

    cudaFree(d_euler);  cudaFree(d_eulerF);
    cudaFree(d_eLcS);   cudaFree(d_eLcF);
    cudaFree(d_eRcS);   cudaFree(d_eRcF);
    cudaFree(d_argmin);
    cudaFree(d_focc);   cudaFree(d_lcount);
    cudaFree(d_tile_num); cudaFree(d_tile_den);
    cudaFreeHost(h_tile_num); cudaFreeHost(h_tile_den);

    env->ReleaseShortArrayElements (j_euler_depths,        h_euler,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(j_euler_F,             h_eulerF, JNI_ABORT);
    env->ReleaseShortArrayElements (j_euler_left_child_s,  h_eLcS,   JNI_ABORT);
    env->ReleaseDoubleArrayElements(j_euler_left_child_f,  h_eLcF,   JNI_ABORT);
    env->ReleaseShortArrayElements (j_euler_right_child_s, h_eRcS,   JNI_ABORT);
    env->ReleaseDoubleArrayElements(j_euler_right_child_f, h_eRcF,   JNI_ABORT);
    env->ReleaseCharArrayElements  (j_sparse_argmin,       h_argmin, JNI_ABORT);
    env->ReleaseIntArrayElements   (j_first_occ,           h_focc,   JNI_ABORT);
    env->ReleaseIntArrayElements   (j_euler_len,           h_elen,   JNI_ABORT);
    env->ReleaseIntArrayElements   (j_leaf_count,          h_lcount, JNI_ABORT);
    output.release();
}

extern "C" JNIEXPORT void JNICALL
Java_stelarx_gpu_GPUSimilarityMatrix_computeSimilarityGPUWide(
    JNIEnv* env,
    jclass cls,
    jintArray    j_euler_depths,
    jdoubleArray j_euler_F,
    jintArray    j_euler_left_child_s,
    jdoubleArray j_euler_left_child_f,
    jintArray    j_euler_right_child_s,
    jdoubleArray j_euler_right_child_f,
    jbyteArray   j_micro_argmin,
    jintArray    j_macro_argmin,
    jintArray    j_first_occ,
    jintArray    j_euler_len,
    jintArray    j_leaf_count,
    jint numTrees,
    jint n,
    jint E_max,
    jint microLog,
    jint blockSize,
    jint blockMax,
    jint macroLog,
    jint tileSizeB,
    jint treeVramCapMiB,
    jdouble progressInterval,
    jint progressMaxSteps,
    jdoubleArray j_num_sum_out,
    jdoubleArray j_den_sum_out,
    jobjectArray j_packed_num_out,
    jobjectArray j_packed_den_out,
    jint packed_segment_shift
) {
    (void)cls;
    jboolean isCopy;
    jint*    h_euler  = env->GetIntArrayElements   (j_euler_depths,         &isCopy);
    jdouble* h_eulerF = env->GetDoubleArrayElements(j_euler_F,              &isCopy);
    jint*    h_eLcS   = env->GetIntArrayElements   (j_euler_left_child_s,   &isCopy);
    jdouble* h_eLcF   = env->GetDoubleArrayElements(j_euler_left_child_f,   &isCopy);
    jint*    h_eRcS   = env->GetIntArrayElements   (j_euler_right_child_s,  &isCopy);
    jdouble* h_eRcF   = env->GetDoubleArrayElements(j_euler_right_child_f,  &isCopy);
    jbyte*   h_micro  = env->GetByteArrayElements  (j_micro_argmin,         &isCopy);
    jint*    h_macro  = env->GetIntArrayElements   (j_macro_argmin,         &isCopy);
    jint*    h_focc   = env->GetIntArrayElements   (j_first_occ,            &isCopy);
    jint*    h_elen   = env->GetIntArrayElements   (j_euler_len,            &isCopy);
    jint*    h_lcount = env->GetIntArrayElements   (j_leaf_count,           &isCopy);
    SimHostOutput output;
    if (!output.acquire(env, j_num_sum_out, j_den_sum_out,
                        j_packed_num_out, j_packed_den_out, packed_segment_shift)) return;

    size_t free_vram = 0, total_vram = 0;
    SIM_CUDA_CHECK(cudaMemGetInfo(&free_vram, &total_vram));

    int B = tileSizeB;
    if (B <= 0) {
        B = (int)ceil(sqrt((double)n * numTrees));
        if (B > n) B = n;
    }
    while (B > 1 && 2LL * B * B * (long long)sizeof(double)
            > (long long)(free_vram * 0.40)) B /= 2;
    if (B < 1) B = 1;

    size_t tile_vram = 2ULL * B * B * sizeof(double);
    size_t headroom = 64ULL * 1024 * 1024;
    size_t available = (free_vram > tile_vram + headroom)
                     ? free_vram - tile_vram - headroom
                     : 16ULL * 1024 * 1024;
    size_t requested = (size_t)treeVramCapMiB * 1024 * 1024;
    size_t remaining = (requested < available) ? requested : available;

    size_t per_tree = (size_t)E_max * 36
                    + (size_t)microLog * E_max * sizeof(unsigned char)
                    + (size_t)macroLog * blockMax * sizeof(int)
                    + (size_t)n * sizeof(int)
                    + sizeof(int);
    int delta = (per_tree > 0) ? (int)(remaining / per_tree) : numTrees;
    if (delta < 1) delta = 1;
    if (delta > numTrees) delta = numTrees;

    int num_batches = (numTrees + delta - 1) / delta;
    int num_tiles_side = (n + B - 1) / B;
    int num_tiles = num_tiles_side * (num_tiles_side + 1) / 2;
    fprintf(stderr,
        "\n[STELAR-X sim] GPU similarity matrix (wide blocked RMQ): n=%d  k=%d  "
        "tile B=%d  tree-batch Δ=%d  (%d batches × %d tiles)\n",
        n, numTrees, B, delta, num_batches, num_tiles);
    fprintf(stderr,
        "[STELAR-X sim] GPU VRAM: tile %.1f MB  tree-data %.1f MB  "
        "(cap %d MiB; free %.0f MB / total %.0f MB)\n",
        tile_vram / 1e6, (double)delta * per_tree / 1e6,
        treeVramCapMiB, free_vram / 1e6, total_vram / 1e6);

    int* d_euler = nullptr;
    double* d_eulerF = nullptr;
    int* d_eLcS = nullptr;
    double* d_eLcF = nullptr;
    int* d_eRcS = nullptr;
    double* d_eRcF = nullptr;
    unsigned char* d_micro = nullptr;
    int* d_macro = nullptr;
    int* d_focc = nullptr;
    int* d_lcount = nullptr;

    long long sz_int_e = (long long)delta * E_max * sizeof(int);
    long long sz_double_e = (long long)delta * E_max * sizeof(double);
    long long sz_micro = (long long)delta * microLog * E_max;
    long long sz_macro = (long long)delta * macroLog * blockMax * sizeof(int);
    SIM_CUDA_CHECK(cudaMalloc(&d_euler, sz_int_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eulerF, sz_double_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eLcS, sz_int_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eLcF, sz_double_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eRcS, sz_int_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_eRcF, sz_double_e));
    SIM_CUDA_CHECK(cudaMalloc(&d_micro, sz_micro));
    SIM_CUDA_CHECK(cudaMalloc(&d_macro, sz_macro));
    SIM_CUDA_CHECK(cudaMalloc(&d_focc, (long long)delta * n * sizeof(int)));
    SIM_CUDA_CHECK(cudaMalloc(&d_lcount, delta * sizeof(int)));

    double* d_tile_num = nullptr;
    double* d_tile_den = nullptr;
    SIM_CUDA_CHECK(cudaMalloc(&d_tile_num, (long long)B * B * sizeof(double)));
    SIM_CUDA_CHECK(cudaMalloc(&d_tile_den, (long long)B * B * sizeof(double)));
    double* h_tile_num = nullptr;
    double* h_tile_den = nullptr;
    SIM_CUDA_CHECK(cudaMallocHost(&h_tile_num, (long long)B * B * sizeof(double)));
    SIM_CUDA_CHECK(cudaMallocHost(&h_tile_den, (long long)B * B * sizeof(double)));

    int total_work = num_batches * num_tiles;
    int work_done = 0;
    double t_start = sim_now_sec();
    double t_last_print = t_start - progressInterval;
    double last_pct_printed = -1.0;
    const bool step_mode = (progressMaxSteps > 0);
    int use_color = sim_use_color();

    for (int t0 = 0; t0 < numTrees; t0 += delta) {
        int dt = (t0 + delta > numTrees) ? numTrees - t0 : delta;
        long long off_e = (long long)t0 * E_max;
        long long bytes_int_e = (long long)dt * E_max * sizeof(int);
        long long bytes_double_e = (long long)dt * E_max * sizeof(double);
        long long off_micro = (long long)t0 * microLog * E_max;
        long long bytes_micro = (long long)dt * microLog * E_max;
        long long off_macro = (long long)t0 * macroLog * blockMax;
        long long bytes_macro = (long long)dt * macroLog * blockMax * sizeof(int);

        SIM_CUDA_CHECK(cudaMemcpy(d_euler, h_euler + off_e,
            bytes_int_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eulerF, h_eulerF + off_e,
            bytes_double_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eLcS, h_eLcS + off_e,
            bytes_int_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eLcF, h_eLcF + off_e,
            bytes_double_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eRcS, h_eRcS + off_e,
            bytes_int_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_eRcF, h_eRcF + off_e,
            bytes_double_e, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_micro,
            reinterpret_cast<unsigned char*>(h_micro) + off_micro,
            bytes_micro, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_macro, h_macro + off_macro,
            bytes_macro, cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_focc, h_focc + (long long)t0 * n,
            (long long)dt * n * sizeof(int), cudaMemcpyHostToDevice));
        SIM_CUDA_CHECK(cudaMemcpy(d_lcount, h_lcount + t0,
            dt * sizeof(int), cudaMemcpyHostToDevice));

        for (int a0 = 0; a0 < n; a0 += B) {
            int bA = (a0 + B > n) ? n - a0 : B;
            for (int b0 = a0; b0 < n; b0 += B) {
                int bB_tile = (b0 + B > n) ? n - b0 : B;
                SIM_CUDA_CHECK(cudaMemset(d_tile_num, 0,
                    (long long)bA * bB_tile * sizeof(double)));
                SIM_CUDA_CHECK(cudaMemset(d_tile_den, 0,
                    (long long)bA * bB_tile * sizeof(double)));

                dim3 block(32, 32);
                dim3 grid((bA + 31) / 32, (bB_tile + 31) / 32);
                sim_tile_kernel_wide<<<grid, block>>>(
                    d_euler, d_eulerF,
                    d_eLcS, d_eLcF, d_eRcS, d_eRcF,
                    d_micro, d_macro, d_focc, d_lcount,
                    dt, n, E_max, microLog, blockSize, blockMax, macroLog,
                    a0, b0, bA, bB_tile, d_tile_num, d_tile_den);
                SIM_CUDA_CHECK(cudaDeviceSynchronize());
                SIM_CUDA_CHECK(cudaMemcpy(h_tile_num, d_tile_num,
                    (long long)bA * bB_tile * sizeof(double), cudaMemcpyDeviceToHost));
                SIM_CUDA_CHECK(cudaMemcpy(h_tile_den, d_tile_den,
                    (long long)bA * bB_tile * sizeof(double), cudaMemcpyDeviceToHost));

                if (!output.packed) {
                    for (int da = 0; da < bA; da++) {
                        int a = a0 + da;
                        for (int db = 0; db < bB_tile; db++) {
                            int b = b0 + db;
                            if (a >= b) continue;
                            long tidx = (long)da * bB_tile + db;
                            double nv = h_tile_num[tidx];
                            double dv = h_tile_den[tidx];
                            if (dv == 0.0) continue;
                            output.dense_num[a * n + b] += nv;
                            output.dense_num[b * n + a] += nv;
                            output.dense_den[a * n + b] += dv;
                            output.dense_den[b * n + a] += dv;
                        }
                    }
                } else {
                    for (int da = 0; da < bA; da++) {
                        int a = a0 + da;
                        for (int db = 0; db < bB_tile; db++) {
                            int b = b0 + db;
                            if (a >= b) continue;
                            long tidx = (long)da * bB_tile + db;
                            double nv = h_tile_num[tidx];
                            double dv = h_tile_den[tidx];
                            if (dv == 0.0) continue;
                            output.addPacked(a, b, n, nv, dv);
                        }
                    }
                }

                work_done++;
                double now = sim_now_sec();
                double pct = (total_work > 0) ? 100.0 * work_done / total_work : 100.0;
                bool is_last = (work_done == total_work);
                bool should_print = step_mode
                    ? (is_last || pct - last_pct_printed >= 100.0 / progressMaxSteps)
                    : (is_last || now - t_last_print >= progressInterval);
                if (should_print) {
                    sim_print_progress(work_done, total_work, now - t_start, use_color, is_last);
                    t_last_print = now;
                    last_pct_printed = pct;
                }
            }
        }
    }

    cudaFree(d_euler); cudaFree(d_eulerF);
    cudaFree(d_eLcS); cudaFree(d_eLcF);
    cudaFree(d_eRcS); cudaFree(d_eRcF);
    cudaFree(d_micro); cudaFree(d_macro);
    cudaFree(d_focc); cudaFree(d_lcount);
    cudaFree(d_tile_num); cudaFree(d_tile_den);
    cudaFreeHost(h_tile_num); cudaFreeHost(h_tile_den);

    env->ReleaseIntArrayElements   (j_euler_depths,        h_euler,  JNI_ABORT);
    env->ReleaseDoubleArrayElements(j_euler_F,             h_eulerF, JNI_ABORT);
    env->ReleaseIntArrayElements   (j_euler_left_child_s,  h_eLcS,   JNI_ABORT);
    env->ReleaseDoubleArrayElements(j_euler_left_child_f,  h_eLcF,   JNI_ABORT);
    env->ReleaseIntArrayElements   (j_euler_right_child_s, h_eRcS,   JNI_ABORT);
    env->ReleaseDoubleArrayElements(j_euler_right_child_f, h_eRcF,   JNI_ABORT);
    env->ReleaseByteArrayElements  (j_micro_argmin,        h_micro,  JNI_ABORT);
    env->ReleaseIntArrayElements   (j_macro_argmin,        h_macro,  JNI_ABORT);
    env->ReleaseIntArrayElements   (j_first_occ,           h_focc,   JNI_ABORT);
    env->ReleaseIntArrayElements   (j_euler_len,           h_elen,   JNI_ABORT);
    env->ReleaseIntArrayElements   (j_leaf_count,          h_lcount, JNI_ABORT);
    output.release();
}
