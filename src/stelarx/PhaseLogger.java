package stelarx;

import stelarx.gpu.GPUWeightCalculator;
import java.util.concurrent.atomic.AtomicLong;

/**
 * Step-by-step phase logging with both per-phase and run-wide peak CPU RAM and
 * GPU VRAM.
 *
 * Two background threads run for every phase:
 *   - cpu-ram-poller : reads /proc/self/status VmRSS every 50 ms (total process RSS,
 *                      includes JNI native memory — matches what htop shows)
 *   - vram-poller    : polls cudaMemGetInfo() every 50 ms (GPU phases only)
 *
 * Output format:
 *   ▶  Phase N  Name  [GPU/CPU]
 *   ...kernel/Java output interleaved...
 *        ✓  1234 ms
 *           peak RAM (this step): 8192 MiB │ peak RAM (so far): 8192 MiB
 *           peak VRAM (this step): 231 MiB │ peak VRAM (so far): 231 MiB
 */
public class PhaseLogger {
    private static volatile String currentPhase = "startup";

    // ── ANSI colours (shared detection with Banner) ───────────────────────────
    private static final boolean COLOR = Banner.useColor();
    private static final String RST  = "\033[0m";
    private static final String BOLD = "\033[1m";
    private static final String DIM  = "\033[2m";
    private static final String CYAN = "\033[36m";
    private static final String GRN  = "\033[32m";
    private static final String YLW  = "\033[33m";
    private static final String WHT  = "\033[97m";

    private static String c(String code, String text) {
        return COLOR ? code + text + RST : text;
    }

    // ── CPU RAM poller (/proc/self/status VmRSS) ──────────────────────────────
    private static final AtomicLong peakCpuRssMiB      = new AtomicLong(0);
    private static final AtomicLong peakCpuRssSoFarMiB = new AtomicLong(0);
    private static volatile boolean cpuPolling       = false;
    private static Thread           cpuPollThread    = null;

    /** Read resident set size in MiB from /proc/self/status, or -1 on failure. */
    private static long readVmRssMiB() {
        try (java.io.BufferedReader br = new java.io.BufferedReader(
                new java.io.FileReader("/proc/self/status"))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("VmRSS:")) {
                    // format: "VmRSS:   12345 kB"
                    String[] parts = line.trim().split("\\s+");
                    if (parts.length >= 2) return Long.parseLong(parts[1]) / 1024;
                }
            }
        } catch (Exception ignored) {}
        // fallback: Java heap only
        Runtime rt = Runtime.getRuntime();
        return (rt.totalMemory() - rt.freeMemory()) / (1024 * 1024);
    }

    private static void startCpuPolling() {
        peakCpuRssMiB.set(0);
        cpuPolling = true;
        cpuPollThread = new Thread(() -> {
            while (cpuPolling) {
                long rss = readVmRssMiB();
                if (rss > 0) {
                    long cur;
                    do { cur = peakCpuRssMiB.get(); }
                    while (rss > cur && !peakCpuRssMiB.compareAndSet(cur, rss));
                }
                try { Thread.sleep(50); } catch (InterruptedException e) { break; }
            }
        }, "cpu-ram-poller");
        cpuPollThread.setDaemon(true);
        cpuPollThread.start();
    }

    /** Stop CPU polling and format the current-step and run-wide RSS peaks. */
    private static String stopCpuPolling() {
        if (cpuPollThread == null) return null;
        cpuPolling = false;
        try { cpuPollThread.join(500); } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        cpuPollThread = null;

        long peak = peakCpuRssMiB.get();
        if (peak <= 0) return null;
        long peakSoFar = updateMaximum(peakCpuRssSoFarMiB, peak);
        // colour: cyan (CPU colour matches [CPU] tag)
        return c(DIM, "peak ") + c(CYAN, "RAM") + c(DIM, " (this step): ")
            + c(CYAN, peak + " MiB")
            + "  " + c(DIM, "│") + "  "
            + c(DIM, "peak ") + c(CYAN, "RAM") + c(DIM, " (so far): ")
            + c(CYAN, peakSoFar + " MiB");
    }

    // ── Background VRAM poller ────────────────────────────────────────────────
    private static final AtomicLong minFreeVRAM = new AtomicLong(Long.MAX_VALUE);
    private static final AtomicLong totalVRAM   = new AtomicLong(0);
    private static final AtomicLong peakVramUsedSoFarMiB = new AtomicLong(0);
    private static volatile boolean vramPolling = false;
    private static Thread           vramPollThread = null;

    private static void startVramPolling() {
        if (!GPUWeightCalculator.isLoaded()) return;
        minFreeVRAM.set(Long.MAX_VALUE);
        totalVRAM.set(0);
        vramPolling = true;
        vramPollThread = new Thread(() -> {
            while (vramPolling) {
                long[] v = GPUWeightCalculator.queryVRAMMiB();
                if (v != null) {
                    long free  = v[0];
                    long total = v[1];
                    totalVRAM.set(total);
                    long cur;
                    do { cur = minFreeVRAM.get(); }
                    while (free < cur && !minFreeVRAM.compareAndSet(cur, free));
                }
                try { Thread.sleep(50); } catch (InterruptedException e) { break; }
            }
        }, "vram-poller");
        vramPollThread.setDaemon(true);
        vramPollThread.start();
    }

    /** Stop VRAM polling and format the current-step and run-wide usage peaks. */
    private static String stopVramPolling() {
        if (vramPollThread == null) return null;
        vramPolling = false;
        try { vramPollThread.join(500); } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        vramPollThread = null;

        long minFree = minFreeVRAM.get();
        long total   = totalVRAM.get();
        if (minFree == Long.MAX_VALUE || total == 0) return null;

        long peakUsed = total - minFree;
        long peakSoFar = updateMaximum(peakVramUsedSoFarMiB, peakUsed);
        String stepCol = (peakUsed * 2 < total) ? GRN : YLW;
        String runCol  = (peakSoFar * 2 < total) ? GRN : YLW;
        return c(DIM, "peak ") + c(GRN, "VRAM") + c(DIM, " (this step): ")
            + c(stepCol, peakUsed + " MiB")
            + "  " + c(DIM, "│") + "  "
            + c(DIM, "peak ") + c(GRN, "VRAM") + c(DIM, " (so far): ")
            + c(runCol, peakSoFar + " MiB");
    }

    /** Atomically raise a run-wide high-water mark and return its new value. */
    private static long updateMaximum(AtomicLong maximum, long candidate) {
        long current;
        do {
            current = maximum.get();
            if (candidate <= current) return current;
        } while (!maximum.compareAndSet(current, candidate));
        return candidate;
    }

    // ── Public API ────────────────────────────────────────────────────────────

    /**
     * Print phase-start header and return a start timestamp.
     * @param label human-readable phase label
     * @param gpu   true if this phase executes on the GPU
     * @return System.nanoTime() snapshot for passing to {@link #end}
     */
    public static long begin(String label, boolean gpu) {
        currentPhase = label;
        startCpuPolling();
        if (gpu) startVramPolling();
        String tag = gpu ? c(GRN, "[GPU]") : c(CYAN, "[CPU]");
        System.err.println();
        System.err.println("  " + c(BOLD, "▶  " + label) + "  " + tag);
        return System.nanoTime();
    }

    /**
     * Print a phase-completion block with time, RAM peaks, and (for GPU) VRAM peaks.
     * @param label same label passed to {@link #begin}
     * @param t0    timestamp returned by {@link #begin}
     * @param gpu   true if this phase ran on the GPU
     */
    public static void end(String label, long t0, boolean gpu) {
        long ms = (System.nanoTime() - t0) / 1_000_000;
        String cpuStr  = stopCpuPolling();
        String vramStr = gpu ? stopVramPolling() : null;

        System.err.println("     " + c(DIM, "✓") + "  " + c(YLW, ms + " ms"));
        if (cpuStr != null) {
            System.err.println("        " + cpuStr);
        }
        if (vramStr != null) {
            System.err.println("        " + vramStr);
        }
        currentPhase = "between phases";
    }

    /** Print the four-field, wrapper-independent summary for a successful analysis. */
    public static void printRunSummary(String tripletScore, long elapsedMs,
                                       boolean gpuExecution) {
        long cpuMiB = readVmRssMiB();
        if (cpuMiB > 0) updateMaximum(peakCpuRssSoFarMiB, cpuMiB);
        long peakCpuMiB = peakCpuRssSoFarMiB.get();
        long peakGpuMiB = peakVramUsedSoFarMiB.get();

        String cpuValue = peakCpuMiB > 0
            ? c(CYAN, peakCpuMiB + " MiB")
            : c(DIM, "N/A");
        String gpuValue = peakGpuMiB > 0
            ? c(GRN, peakGpuMiB + " MiB")
            : c(DIM, gpuExecution ? "N/A (not observed)" : "N/A (CPU execution)");

        System.err.println();
        System.err.println("  " + c(BOLD, "Run Summary") + "  " + c(DIM, "─".repeat(50)));
        System.err.println();
        System.err.println("    " + summaryRow("Triplet score", c(WHT, tripletScore)));
        System.err.println("    " + summaryRow("Running time", c(YLW,
            String.format(java.util.Locale.ROOT, "%.3f s", elapsedMs / 1000.0))));
        System.err.println("    " + summaryRow("Max CPU RAM", cpuValue));
        System.err.println("    " + summaryRow("Max GPU VRAM", gpuValue));
        System.err.println();
    }

    private static String summaryRow(String label, String value) {
        String paddedLabel = String.format(java.util.Locale.ROOT, "%-18s", label);
        return c(DIM, paddedLabel) + " " + value;
    }

    /** Last phase label, included in fatal diagnostic reports. */
    public static String currentPhase() { return currentPhase; }
}
