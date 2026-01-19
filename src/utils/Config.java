package utils;

/**
 * Configuration class for wQFM-TREE algorithm settings.
 * This class holds global configuration parameters and enums used throughout
 * the algorithm.
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
     * Computation mode for weight calculation (internal use)
     */
    public enum ComputationMode {
        CPU_SINGLE, // Single-threaded CPU computation (memory-optimized)
        CPU_PARALLEL, // Multi-threaded CPU computation (memory-optimized)
        GPU_PARALLEL // GPU-accelerated computation (memory-optimized)
    }

    /**
     * Current computation mode setting (derived from USE_GPU and NUM_THREADS)
     */
    public static ComputationMode COMPUTATION_MODE = ComputationMode.CPU_PARALLEL;

    /**
     * Whether to use GPU acceleration
     * Default: true (GPU is default, use --cpu to disable)
     */
    public static boolean USE_GPU = true;

    /**
     * Number of threads to use for CPU parallelism (-t flag)
     * Default: 0 means use all available processors
     */
    public static int NUM_THREADS = 0;

    /**
     * Get the effective number of threads to use.
     * 
     * @return Number of threads (all available if NUM_THREADS is 0)
     */
    public static int getEffectiveThreadCount() {
        if (NUM_THREADS <= 0) {
            return Runtime.getRuntime().availableProcessors();
        }
        return NUM_THREADS;
    }

    /**
     * Update the computation mode based on USE_GPU and NUM_THREADS settings.
     * Should be called after parsing command line arguments.
     */
    public static void updateComputationMode() {
        if (USE_GPU) {
            COMPUTATION_MODE = ComputationMode.GPU_PARALLEL;
        } else if (getEffectiveThreadCount() == 1) {
            COMPUTATION_MODE = ComputationMode.CPU_SINGLE;
        } else {
            COMPUTATION_MODE = ComputationMode.CPU_PARALLEL;
        }
    }
}
