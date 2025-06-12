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
}
