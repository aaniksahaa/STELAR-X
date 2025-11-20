package utils;

/**
 * This is under-development
 */
public class BipartitionExpansionConfig {
    
    /**
     * Bipartition expansion methods
     */
    public enum ExpansionMethod {
        NONE,              // No expansion - only gene tree bipartitions
        DISTANCE_ONLY,     // Only distance-based expansion
        CONSENSUS_ONLY,    // Only consensus-based expansion  
        DISTANCE_CONSENSUS,// Both distance and consensus
        FULL              // All methods including polytomy resolution
    }
    
    /**
     * Distance matrix methods for tree inference
     */
    public enum DistanceMethod {
        UPGMA,            // UPGMA clustering
        NEIGHBOR_JOINING, // Neighbor-joining
        BOTH             // Try both methods
    }
    
    /**
     * Current expansion method setting
     */
    public static ExpansionMethod EXPANSION_METHOD = ExpansionMethod.DISTANCE_CONSENSUS;
    
    /**
     * Distance method for distance-based expansion
     */
    public static DistanceMethod DISTANCE_METHOD = DistanceMethod.UPGMA;
    
    /**
     * Enable/disable distance-based bipartition expansion
     */
    public static boolean ENABLE_DISTANCE_EXPANSION = true;
    
    /**
     * Enable/disable consensus-based bipartition expansion
     */
    public static boolean ENABLE_CONSENSUS_EXPANSION = true;
    
    /**
     * Enable/disable polytomy resolution
     */
    public static boolean ENABLE_POLYTOMY_RESOLUTION = false;
    
    /**
     * Consensus support thresholds
     * Lower values are more permissive, higher values are more strict
     */
    public static double[] CONSENSUS_THRESHOLDS = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
    
    /**
     * Maximum number of consensus trees to generate per threshold
     */
    public static int MAX_CONSENSUS_TREES = 10;
    
    /**
     * Maximum polytomy size to resolve (for computational efficiency)
     * Polytomies larger than this will be skipped
     */
    public static int MAX_POLYTOMY_SIZE = 10;
    
    /**
     * Number of sampling rounds for polytomy resolution
     */
    public static int POLYTOMY_SAMPLING_ROUNDS = 5;
    
    /**
     * Minimum frequency threshold for adding bipartitions from polytomy resolution
     */
    public static int MIN_BIPARTITION_FREQUENCY = 2;
    
    /**
     * Enable verbose output for expansion process
     */
    public static boolean VERBOSE_EXPANSION = false;
    
    /**
     * Maximum number of additional bipartitions to add (safety limit)
     * Set to -1 to disable limit (no trimming)
     */
    public static int MAX_ADDITIONAL_BIPARTITIONS = -1;
    
    /**
     * Check if any expansion method is enabled
     */
    public static boolean isExpansionEnabled() {
        return EXPANSION_METHOD != ExpansionMethod.NONE;
    }
    
    /**
     * Check if distance-based expansion is enabled
     */
    public static boolean isDistanceExpansionEnabled() {
        return ENABLE_DISTANCE_EXPANSION && 
               (EXPANSION_METHOD == ExpansionMethod.DISTANCE_ONLY ||
                EXPANSION_METHOD == ExpansionMethod.DISTANCE_CONSENSUS ||
                EXPANSION_METHOD == ExpansionMethod.FULL);
    }
    
    /**
     * Check if consensus-based expansion is enabled
     */
    public static boolean isConsensusExpansionEnabled() {
        return ENABLE_CONSENSUS_EXPANSION &&
               (EXPANSION_METHOD == ExpansionMethod.CONSENSUS_ONLY ||
                EXPANSION_METHOD == ExpansionMethod.DISTANCE_CONSENSUS ||
                EXPANSION_METHOD == ExpansionMethod.FULL);
    }
    
    /**
     * Check if polytomy resolution is enabled
     */
    public static boolean isPolytomyResolutionEnabled() {
        return ENABLE_POLYTOMY_RESOLUTION &&
               EXPANSION_METHOD == ExpansionMethod.FULL;
    }
    
    /**
     * Print current configuration
     */
    public static void printConfig() {
        System.out.println("Bipartition Expansion Configuration:");
        System.out.println("  Expansion Method: " + EXPANSION_METHOD);
        System.out.println("  Distance Method: " + DISTANCE_METHOD);
        System.out.println("  Distance Expansion: " + isDistanceExpansionEnabled());
        System.out.println("  Consensus Expansion: " + isConsensusExpansionEnabled());
        System.out.println("  Polytomy Resolution: " + isPolytomyResolutionEnabled());
        
        if (isConsensusExpansionEnabled()) {
            System.out.print("  Consensus Thresholds: [");
            for (int i = 0; i < CONSENSUS_THRESHOLDS.length; i++) {
                if (i > 0) System.out.print(", ");
                System.out.print(CONSENSUS_THRESHOLDS[i]);
            }
            System.out.println("]");
            System.out.println("  Max Consensus Trees: " + MAX_CONSENSUS_TREES);
        }
        
        if (isPolytomyResolutionEnabled()) {
            System.out.println("  Max Polytomy Size: " + MAX_POLYTOMY_SIZE);
            System.out.println("  Polytomy Sampling Rounds: " + POLYTOMY_SAMPLING_ROUNDS);
        }
        
        if (MAX_ADDITIONAL_BIPARTITIONS == -1) {
            System.out.println("  Max Additional Bipartitions: unlimited");
        } else {
            System.out.println("  Max Additional Bipartitions: " + MAX_ADDITIONAL_BIPARTITIONS);
        }
        System.out.println("  Verbose Output: " + VERBOSE_EXPANSION);
    }
}
