package stelarx.completion;

import stelarx.util.Threading;

import java.util.Arrays;

/**
 * Builds sorted-row index arrays from a flat n×n distance matrix.
 *
 * Output: int[] sortedRows of size n*n, where
 *   sortedRows[x * n + rank] = taxon ID of x's rank-th nearest neighbor
 *   (rank 0 = x itself at distance 0; rank 1 = nearest other taxon).
 *
 * Used by TreeCompleter's four-point completion algorithm to find, for each
 * missing taxon x, the closest currently-present taxon (the "anchor").
 *
 * Reference: follows the sortByDistance / assureOrderedTaxa ordering of
 * ASTRAL-MP's AbstractMatrix so completed trees are comparable.
 *
 * HINT for n-ary extension: no changes needed here — sortedRows is purely a
 * property of the distance matrix and has nothing to do with tree topology.
 */
public class SortedRowsBuilder {

    /**
     * CPU path: sort every row of the distance matrix in parallel.
     *
     * Time:   O(n² log n)
     * Memory: O(n²) for the output array + O(n) scratch per thread
     *
     * @param dist flat n×n distance matrix, row-major (dist[a*n+b] = D(a,b))
     * @param n    number of taxa
     * @return int[] of size n*n; sortedRows[x*n + rank] = taxon ID at that rank
     */
    public static int[] buildCPU(double[] dist, int n) {
        int[] sortedRows = new int[n * n];

        // Each row is independent: safe to parallelise across rows.
        Threading.processRangeParallel(n, x -> {
            // Build an index array [0, 1, ..., n-1] for this row.
            Integer[] indices = new Integer[n];
            for (int j = 0; j < n; j++) indices[j] = j;

            final int base = x * n;
            // Sort ascending by distance (= descending similarity).
            // Tie-break: descending taxon ID, matching ASTRAL-MP's AbstractMatrix.sortColumn
            // which does: comp==0 ? -o1.compareTo(o2) : comp  (higher ID wins on tie).
            Arrays.sort(indices, (a, b) -> {
                int c = Double.compare(dist[base + a], dist[base + b]);
                return c != 0 ? c : Integer.compare(b, a);  // descending ID on tie
            });

            // Write into flat output array.
            for (int rank = 0; rank < n; rank++) {
                sortedRows[base + rank] = indices[rank];
            }
        });

        return sortedRows;
    }

    /**
     * Exact segmented-row variant used when a flat n*n int array is impossible.
     * Comparator and tie order are identical to {@link #buildCPU(double[], int)}.
     */
    public static int[][] buildPackedCPU(SimilarityMatrix sim) {
        int n = sim.n;
        int[][] sortedRows = new int[n][];
        Threading.processRangeParallel(n, x -> {
            Integer[] indices = new Integer[n];
            for (int j = 0; j < n; j++) indices[j] = j;
            Arrays.sort(indices, (a, b) -> {
                int c = Double.compare(sim.getDist(x, a), sim.getDist(x, b));
                return c != 0 ? c : Integer.compare(b, a);
            });
            int[] row = new int[n];
            for (int rank = 0; rank < n; rank++) row[rank] = indices[rank];
            sortedRows[x] = row;
        });
        return sortedRows;
    }

    // TODO (GPU path): implement buildGPU(double[] dist, int n) using batched
    // Thrust argsort as described in the implementation plan Section 3.2.
    // Signature: public static int[] buildGPU(double[] dist, int n)
    // Native method: native void computeSortedRowsGPU(double[] dist, int n, int batchSize, int[] sortedRows)
}
