package utils;

/**
 * Configuration class for wQFM-TREE algorithm settings.
 * This class holds global configuration parameters and enums used throughout the algorithm.
 */
public class Config {
    
    /**
     * Number of spaces for output formatting (dummy config for initial setup)
     */
    public static int n = 4;
    
    /**
     * Whether to resolve polytomies in gene trees
     */
    public static boolean RESOLVE_POLYTOMY = false;
    
    /**
     * Non-quartet handling strategy enumeration
     */
    public enum NonQuartetType {
        A, B
    }
    
    /**
     * Current non-quartet type setting
     */
    public static NonQuartetType NON_QUARTET_TYPE = NonQuartetType.A;

    /**
     * Computation mode for weight calculation
     */
    public enum ComputationMode {
        CPU_SINGLE,    // Single-threaded CPU computation (memory-optimized)
        CPU_PARALLEL,  // Multi-threaded CPU computation (memory-optimized)
        GPU_PARALLEL   // GPU-accelerated computation (memory-optimized)
    }

    /**
     * Current computation mode setting
     */
    public static ComputationMode COMPUTATION_MODE = ComputationMode.CPU_PARALLEL;

    /**
     * Scoring mode for bipartition weight calculation.
     * 
     * TRIPLET: Uses triplet matching score (STELAR-X default)
     *          Compares bipartition (A|B) vs gene tree bipartition (X|Y)
     *          Formula: score1 + score2 where score_i = p1 * p2 * (p1 + p2 - 2) / 2
     * 
     * QUARTET: Uses quartet matching score (ASTRAL-style)
     *          Compares tripartition (A|B|C) vs gene tree tripartition (X|Y|Z)
     *          where C = S - A - B (complement) and Z = Sg - X - Y
     *          Uses F(a,b,c) = a*b*c*(a+b+c-3)/2 formula with 6 terms
     */
    public enum ScoringMode {
        TRIPLET,   // STELAR-X triplet matching score (default)
        QUARTET    // ASTRAL-style quartet matching score
    }

    /**
     * Current scoring mode setting (default: TRIPLET for backward compatibility)
     */
    public static ScoringMode SCORING_MODE = ScoringMode.TRIPLET;
}
