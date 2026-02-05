package core;

import tree.RangeBipartition;
import tree.MixedBipartition;

/**
 * Quartet score calculator implementing ASTRAL-style weight computation.
 * 
 * ASTRAL's quartet scoring formula:
 * For tripartition T = (A|B|C) vs gene tree node partition M = (X|Y|Z):
 * 
 * 1. Compute 3x3 intersection grid n[i][j] = |side_i of T ∩ side_j of M|
 * 2. F(a,b,c) = a * b * c * (a + b + c - 3) / 2
 * 3. QI(T,M) = sum of 6 terms: F(n11,n22,n33) + F(n11,n23,n32) + 
 *              F(n12,n21,n33) + F(n12,n23,n31) + F(n13,n21,n32) + F(n13,n22,n31)
 * 
 * Optimization: Only compute 4 intersections (ax, ay, bx, by), derive the rest using:
 * - Row sums: a = |A ∩ Sg|, b = |B ∩ Sg|, c = |Sg| - a - b
 * - Column sums: |X|, |Y|, |Z| = |Sg| - |X| - |Y|
 * - Derivation: az = a - ax - ay, bz = b - bx - by, etc.
 * 
 * C (third side of candidate) is NEVER stored or computed directly - it's implicit.
 */
public class QuartetScoreCalculator {
    
    private final InverseIndexManager inverseIndex;
    private final MissingTaxaManager missingTaxa;
    private final int numTaxa;
    
    // Statistics for performance monitoring
    private long totalQuartetCalculations = 0;
    private long totalIntersectionCalculations = 0;
    private long totalGridReconstructions = 0;
    
    // Debug flag
    private boolean debugEnabled = false;
    
    public QuartetScoreCalculator(InverseIndexManager inverseIndex, MissingTaxaManager missingTaxa, int numTaxa) {
        this.inverseIndex = inverseIndex;
        this.missingTaxa = missingTaxa;
        this.numTaxa = numTaxa;
        
        System.out.println("QuartetScoreCalculator initialized");
        System.out.println("  Using ASTRAL-style F(a,b,c) = a*b*c*(a+b+c-3)/2 formula");
        System.out.println("  Optimized: 4 intersections + row/column sum derivation");
    }
    
    /**
     * Enable debug logging for score calculations.
     */
    public void setDebugEnabled(boolean enabled) {
        this.debugEnabled = enabled;
    }
    
    /**
     * ASTRAL's F function: F(a,b,c) = a * b * c * (a + b + c - 3)
     * 
     * This counts the number of shared induced quartets.
     * Note: No division by 2 - matches ASTRAL's actual implementation.
     */
    private double F(int a, int b, int c) {
        if (a < 0 || b < 0 || c < 0) {
            if (debugEnabled) {
                System.err.println("Warning: Negative value in F(" + a + ", " + b + ", " + c + ")");
            }
            return 0;
        }
        
        long sum = (long) a + b + c;
        if (sum < 3) {
            return 0;
        }
        
        // Use long arithmetic to avoid overflow
        // NO division by 2 - this matches ASTRAL's formula
        return (double) ((long) a * b * c * (sum - 3));
    }
    
    /**
     * Compute the 6-term QI score from a 3x3 intersection grid.
     * 
     * The grid represents:
     *   n[0][0] = |A ∩ X|, n[0][1] = |A ∩ Y|, n[0][2] = |A ∩ Z|
     *   n[1][0] = |B ∩ X|, n[1][1] = |B ∩ Y|, n[1][2] = |B ∩ Z|
     *   n[2][0] = |C ∩ X|, n[2][1] = |C ∩ Y|, n[2][2] = |C ∩ Z|
     * 
     * QI = F(n00,n11,n22) + F(n00,n12,n21) + 
     *      F(n01,n10,n22) + F(n01,n12,n20) +
     *      F(n02,n10,n21) + F(n02,n11,n20)
     */
    private double computeQI(int[][] grid) {
        totalGridReconstructions++;
        
        double qi = F(grid[0][0], grid[1][1], grid[2][2])
                  + F(grid[0][0], grid[1][2], grid[2][1])
                  + F(grid[0][1], grid[1][0], grid[2][2])
                  + F(grid[0][1], grid[1][2], grid[2][0])
                  + F(grid[0][2], grid[1][0], grid[2][1])
                  + F(grid[0][2], grid[1][1], grid[2][0]);
        
        return qi;  // No division - matches ASTRAL's formula
    }
    
    /**
     * Calculate quartet score between a candidate bipartition and a gene tree bipartition.
     * 
     * Candidate bipartition (A|B) induces tripartition (A|B|C) where C = S - A - B.
     * Gene tree bipartition (X|Y) induces tripartition (X|Y|Z) where Z = Sg - X - Y.
     * 
     * We only compute 4 intersections and derive the rest:
     * - ax = |A ∩ X|, ay = |A ∩ Y|
     * - bx = |B ∩ X|, by = |B ∩ Y|
     * - Row sums: a = |A ∩ Sg|, b = |B ∩ Sg|
     * - Column sums: x = |X|, y = |Y|, z = |Sg| - x - y
     * 
     * @param candidate Candidate bipartition (A|B)
     * @param geneTree Gene tree bipartition (X|Y)
     * @return Quartet score contribution (QI/2)
     */
    public double calculateQuartetScore(RangeBipartition candidate, RangeBipartition geneTree) {
        totalQuartetCalculations++;
        
        int geneTreeIdx = geneTree.geneTreeIndex;
        int sgSize = missingTaxa.getPresentCount(geneTreeIdx);
        
        // Get sizes of gene tree bipartition sides
        int xSize = geneTree.leftSize();
        int ySize = geneTree.rightSize();
        int zSize = sgSize - xSize - ySize;
        
        if (zSize < 0) {
            if (debugEnabled) {
                System.err.println("Warning: Negative Z size in gene tree " + geneTreeIdx);
            }
            zSize = 0;
        }
        
        // Get restricted sizes of candidate sides in this gene tree
        // a = |A ∩ Sg|, b = |B ∩ Sg|
        int a = missingTaxa.getRestrictedClusterSize(inverseIndex, 
                candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd, geneTreeIdx);
        int b = missingTaxa.getRestrictedClusterSize(inverseIndex,
                candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd, geneTreeIdx);
        int c = sgSize - a - b;
        
        if (c < 0) {
            if (debugEnabled) {
                System.err.println("Warning: Negative C size for candidate in gene tree " + geneTreeIdx);
            }
            c = 0;
        }
        
        // Compute 4 intersections
        int ax = inverseIndex.getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
            geneTreeIdx, geneTree.leftStart, geneTree.leftEnd);
        
        int ay = inverseIndex.getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.leftStart, candidate.leftEnd,
            geneTreeIdx, geneTree.rightStart, geneTree.rightEnd);
        
        int bx = inverseIndex.getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
            geneTreeIdx, geneTree.leftStart, geneTree.leftEnd);
        
        int by = inverseIndex.getRangeIntersectionSize(
            candidate.geneTreeIndex, candidate.rightStart, candidate.rightEnd,
            geneTreeIdx, geneTree.rightStart, geneTree.rightEnd);
        
        totalIntersectionCalculations += 4;
        
        // Derive the remaining cells using row and column sums
        // Row-wise derivation for A and B with Z
        int az = a - ax - ay;
        int bz = b - bx - by;
        
        // Column-wise derivation for C
        int cx = xSize - ax - bx;
        int cy = ySize - ay - by;
        int cz = zSize - az - bz;
        
        // Sanity checks (in debug mode)
        if (debugEnabled) {
            validateGrid(ax, ay, az, bx, by, bz, cx, cy, cz,
                        a, b, c, xSize, ySize, zSize, sgSize);
        }
        
        // Build grid and compute QI
        int[][] grid = {
            {ax, ay, az},
            {bx, by, bz},
            {cx, cy, cz}
        };
        
        return computeQI(grid);
    }
    
    /**
     * Calculate quartet score for a MixedBipartition (cross-tree recombination).
     * 
     * MixedBipartition allows sides A and B to come from different gene trees.
     * The logic is identical except we use the respective tree indices for each side.
     */
    public double calculateQuartetScore(MixedBipartition mixed, RangeBipartition geneTree) {
        totalQuartetCalculations++;
        
        int geneTreeIdx = geneTree.geneTreeIndex;
        int sgSize = missingTaxa.getPresentCount(geneTreeIdx);
        
        // Get sizes of gene tree bipartition sides
        int xSize = geneTree.leftSize();
        int ySize = geneTree.rightSize();
        int zSize = sgSize - xSize - ySize;
        
        if (zSize < 0) {
            zSize = 0;
        }
        
        // Get restricted sizes - note: use respective tree indices for each side
        int a = missingTaxa.getRestrictedClusterSize(inverseIndex,
                mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd, geneTreeIdx);
        int b = missingTaxa.getRestrictedClusterSize(inverseIndex,
                mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd, geneTreeIdx);
        int c = sgSize - a - b;
        
        if (c < 0) {
            c = 0;
        }
        
        // Compute 4 intersections - use respective tree indices for left/right
        int ax = inverseIndex.getRangeIntersectionSize(
            mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
            geneTreeIdx, geneTree.leftStart, geneTree.leftEnd);
        
        int ay = inverseIndex.getRangeIntersectionSize(
            mixed.leftTreeIndex, mixed.leftStart, mixed.leftEnd,
            geneTreeIdx, geneTree.rightStart, geneTree.rightEnd);
        
        int bx = inverseIndex.getRangeIntersectionSize(
            mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
            geneTreeIdx, geneTree.leftStart, geneTree.leftEnd);
        
        int by = inverseIndex.getRangeIntersectionSize(
            mixed.rightTreeIndex, mixed.rightStart, mixed.rightEnd,
            geneTreeIdx, geneTree.rightStart, geneTree.rightEnd);
        
        totalIntersectionCalculations += 4;
        
        // Derive remaining cells
        int az = a - ax - ay;
        int bz = b - bx - by;
        int cx = xSize - ax - bx;
        int cy = ySize - ay - by;
        int cz = zSize - az - bz;
        
        // Build grid and compute QI
        int[][] grid = {
            {ax, ay, az},
            {bx, by, bz},
            {cx, cy, cz}
        };
        
        return computeQI(grid);
    }
    
    /**
     * Validate that the 3x3 grid is consistent with row and column sums.
     * Used for debugging to catch errors in the derivation.
     */
    private void validateGrid(int ax, int ay, int az, int bx, int by, int bz,
                             int cx, int cy, int cz,
                             int a, int b, int c, int x, int y, int z, int sg) {
        // Check for negative values
        if (ax < 0 || ay < 0 || az < 0 || bx < 0 || by < 0 || bz < 0 ||
            cx < 0 || cy < 0 || cz < 0) {
            System.err.println("Grid validation: NEGATIVE VALUES DETECTED");
            System.err.println("  Grid: [[" + ax + "," + ay + "," + az + "], [" +
                             bx + "," + by + "," + bz + "], [" +
                             cx + "," + cy + "," + cz + "]]");
        }
        
        // Check row sums
        if (ax + ay + az != a) {
            System.err.println("Grid validation: Row A mismatch: " + (ax + ay + az) + " != " + a);
        }
        if (bx + by + bz != b) {
            System.err.println("Grid validation: Row B mismatch: " + (bx + by + bz) + " != " + b);
        }
        if (cx + cy + cz != c) {
            System.err.println("Grid validation: Row C mismatch: " + (cx + cy + cz) + " != " + c);
        }
        
        // Check column sums
        if (ax + bx + cx != x) {
            System.err.println("Grid validation: Column X mismatch: " + (ax + bx + cx) + " != " + x);
        }
        if (ay + by + cy != y) {
            System.err.println("Grid validation: Column Y mismatch: " + (ay + by + cy) + " != " + y);
        }
        if (az + bz + cz != z) {
            System.err.println("Grid validation: Column Z mismatch: " + (az + bz + cz) + " != " + z);
        }
        
        // Check total
        int total = ax + ay + az + bx + by + bz + cx + cy + cz;
        if (total != sg) {
            System.err.println("Grid validation: Total mismatch: " + total + " != " + sg);
        }
    }
    
    /**
     * Get statistics about quartet score calculations.
     */
    public String getStatistics() {
        StringBuilder sb = new StringBuilder();
        sb.append("Quartet Score Calculator Statistics:\n");
        sb.append("  Total quartet calculations: ").append(totalQuartetCalculations).append("\n");
        sb.append("  Total intersections computed: ").append(totalIntersectionCalculations).append("\n");
        sb.append("  Total grid reconstructions: ").append(totalGridReconstructions).append("\n");
        if (totalQuartetCalculations > 0) {
            sb.append("  Intersections per calculation: ").append(
                totalIntersectionCalculations / (double) totalQuartetCalculations).append("\n");
        }
        return sb.toString();
    }
    
    /**
     * Reset statistics counters.
     */
    public void resetStatistics() {
        totalQuartetCalculations = 0;
        totalIntersectionCalculations = 0;
        totalGridReconstructions = 0;
    }
}
