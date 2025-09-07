package expansion;

import java.util.*;
import preprocessing.GeneTrees;
import tree.STBipartition;
import tree.Tree;
import utils.BitSet;
import utils.BipartitionExpansionConfig;

/**
 * Main manager for bipartition expansion functionality.
 * 
 * This class coordinates all bipartition expansion methods and provides
 * a unified interface for expanding the candidate bipartition set beyond
 * just the bipartitions found in input gene trees.
 * 
 * Based on ASTRAL-II expansion strategies:
 * 1. Distance-based expansion (UPGMA/NJ trees from distance matrices)
 * 2. Consensus-based expansion (greedy consensus with multiple thresholds)
 * 3. Polytomy resolution (advanced method for resolving multi-way splits)
 */
public class BipartitionExpansionManager {
    
    private GeneTrees geneTrees;
    private DistanceMatrix distanceMatrix;
    private ConsensusTreeBuilder consensusBuilder;
    
    // Statistics
    private int originalBipartitionCount;
    private int distanceBasedCount;
    private int consensusBasedCount;
    private int polytomyResolutionCount;
    
    public BipartitionExpansionManager(GeneTrees geneTrees) {
        this.geneTrees = geneTrees;
        this.originalBipartitionCount = geneTrees.stBipartitions.size();
        
        // Initialize expansion components based on configuration
        if (BipartitionExpansionConfig.isDistanceExpansionEnabled()) {
            this.distanceMatrix = new DistanceMatrix(geneTrees);
        }
        
        if (BipartitionExpansionConfig.isConsensusExpansionEnabled()) {
            this.consensusBuilder = new ConsensusTreeBuilder(geneTrees);
        }
    }
    
    /**
     * Main method to expand the candidate bipartition set.
     * This method applies all enabled expansion techniques and returns
     * an expanded list of candidate bipartitions.
     */
    public List<STBipartition> expandBipartitions(List<STBipartition> originalCandidates) {
        if (!BipartitionExpansionConfig.isExpansionEnabled()) {
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("Bipartition expansion is disabled. Using original candidates only.");
            }
            return new ArrayList<>(originalCandidates);
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("\n=== BIPARTITION EXPANSION STATISTICS ===");
            System.out.println("Starting bipartition expansion...");
            System.out.println("Original candidate count: " + originalCandidates.size());
            BipartitionExpansionConfig.printConfig();
            System.out.println("==========================================");
        }
        
        // Start with original candidates
        Set<STBipartition> expandedCandidates = new HashSet<>(originalCandidates);
        
        // Apply distance-based expansion
        if (BipartitionExpansionConfig.isDistanceExpansionEnabled()) {
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("\n--- Distance-Based Expansion ---");
            }
            List<STBipartition> distanceBipartitions = performDistanceBasedExpansion();
            int beforeSize = expandedCandidates.size();
            expandedCandidates.addAll(distanceBipartitions);
            distanceBasedCount = expandedCandidates.size() - beforeSize;
            
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("Distance bipartitions found: " + distanceBipartitions.size());
                System.out.println("New unique distance bipartitions added: " + distanceBasedCount);
                System.out.println("Total candidates after distance expansion: " + expandedCandidates.size());
            }
        }
        
        // Apply consensus-based expansion
        if (BipartitionExpansionConfig.isConsensusExpansionEnabled()) {
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("\n--- Consensus-Based Expansion ---");
            }
            List<STBipartition> consensusBipartitions = performConsensusBasedExpansion();
            int beforeSize = expandedCandidates.size();
            expandedCandidates.addAll(consensusBipartitions);
            consensusBasedCount = expandedCandidates.size() - beforeSize;
            
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("Consensus bipartitions found: " + consensusBipartitions.size());
                System.out.println("New unique consensus bipartitions added: " + consensusBasedCount);
                System.out.println("Total candidates after consensus expansion: " + expandedCandidates.size());
            }
        }
        
        // Apply polytomy resolution (if enabled)
        if (BipartitionExpansionConfig.isPolytomyResolutionEnabled()) {
            List<STBipartition> polytomyBipartitions = performPolytomyResolution();
            int beforeSize = expandedCandidates.size();
            expandedCandidates.addAll(polytomyBipartitions);
            polytomyResolutionCount = expandedCandidates.size() - beforeSize;
            
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("Added " + polytomyResolutionCount + " polytomy resolution bipartitions");
            }
        }
        
        // Apply safety limit (if enabled)
        List<STBipartition> result = new ArrayList<>(expandedCandidates);
        if (BipartitionExpansionConfig.MAX_ADDITIONAL_BIPARTITIONS != -1 && 
            result.size() > BipartitionExpansionConfig.MAX_ADDITIONAL_BIPARTITIONS) {
            
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.out.println("Applying safety limit: reducing from " + result.size() + 
                                 " to " + BipartitionExpansionConfig.MAX_ADDITIONAL_BIPARTITIONS);
            }
            
            // Keep original candidates and top additional ones
            result = new ArrayList<>(originalCandidates);
            Set<STBipartition> originalSet = new HashSet<>(originalCandidates);
            
            List<STBipartition> additionalCandidates = new ArrayList<>();
            for (STBipartition bip : expandedCandidates) {
                if (!originalSet.contains(bip)) {
                    additionalCandidates.add(bip);
                }
            }
            
            // Add additional candidates up to the limit
            int remainingSlots = BipartitionExpansionConfig.MAX_ADDITIONAL_BIPARTITIONS - originalCandidates.size();
            for (int i = 0; i < Math.min(remainingSlots, additionalCandidates.size()); i++) {
                result.add(additionalCandidates.get(i));
            }
        } else if (BipartitionExpansionConfig.VERBOSE_EXPANSION && BipartitionExpansionConfig.MAX_ADDITIONAL_BIPARTITIONS == -1) {
            System.out.println("No safety limit applied (unlimited bipartitions allowed)");
        }
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            printExpansionSummary(originalCandidates.size(), result.size());
            System.out.println("=========================================\n");
        }
        
        return result;
    }
    
    /**
     * Perform distance-based bipartition expansion
     */
    private List<STBipartition> performDistanceBasedExpansion() {
        List<STBipartition> distanceBipartitions = new ArrayList<>();
        
        if (distanceMatrix == null) {
            return distanceBipartitions;
        }
        
        try {
            // Build distance-based tree(s)
            if (BipartitionExpansionConfig.DISTANCE_METHOD == BipartitionExpansionConfig.DistanceMethod.UPGMA ||
                BipartitionExpansionConfig.DISTANCE_METHOD == BipartitionExpansionConfig.DistanceMethod.BOTH) {
                
                Tree upgmaTree = distanceMatrix.buildUPGMATree();
                if (upgmaTree != null) {
                    List<STBipartition> upgmaBipartitions = distanceMatrix.extractBipartitions(upgmaTree);
                    distanceBipartitions.addAll(upgmaBipartitions);
                }
            }
            
            // TODO: Add Neighbor-Joining implementation if needed
            if (BipartitionExpansionConfig.DISTANCE_METHOD == BipartitionExpansionConfig.DistanceMethod.NEIGHBOR_JOINING ||
                BipartitionExpansionConfig.DISTANCE_METHOD == BipartitionExpansionConfig.DistanceMethod.BOTH) {
                // Neighbor-joining implementation would go here
                if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                    System.out.println("Neighbor-joining not yet implemented, using UPGMA only.");
                }
            }
            
        } catch (Exception e) {
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.err.println("Error in distance-based expansion: " + e.getMessage());
            }
        }
        
        return distanceBipartitions;
    }
    
    /**
     * Perform consensus-based bipartition expansion
     */
    private List<STBipartition> performConsensusBasedExpansion() {
        List<STBipartition> consensusBipartitions = new ArrayList<>();
        
        if (consensusBuilder == null) {
            return consensusBipartitions;
        }
        
        try {
            // Build consensus trees with multiple thresholds
            List<Tree> consensusTrees = consensusBuilder.buildConsensusTreesWithThresholds();
            
            // Extract bipartitions from all consensus trees
            consensusBipartitions = consensusBuilder.extractBipartitionsFromConsensusTrees(consensusTrees);
            
        } catch (Exception e) {
            if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
                System.err.println("Error in consensus-based expansion: " + e.getMessage());
            }
        }
        
        return consensusBipartitions;
    }
    
    /**
     * Perform polytomy resolution (advanced method)
     */
    private List<STBipartition> performPolytomyResolution() {
        List<STBipartition> polytomyBipartitions = new ArrayList<>();
        
        // This is a placeholder for polytomy resolution
        // The full implementation would be quite complex and similar to ASTRAL's approach
        
        if (BipartitionExpansionConfig.VERBOSE_EXPANSION) {
            System.out.println("Polytomy resolution is not yet fully implemented.");
        }
        
        // TODO: Implement polytomy resolution
        // 1. Identify polytomies in consensus trees
        // 2. Use distance-based methods to resolve them
        // 3. Use sampling-based methods for additional resolutions
        // 4. Extract bipartitions from resolved trees
        
        return polytomyBipartitions;
    }
    
    /**
     * Validate that bipartitions are well-formed
     */
    private List<STBipartition> validateBipartitions(List<STBipartition> bipartitions) {
        List<STBipartition> validBipartitions = new ArrayList<>();
        
        for (STBipartition bip : bipartitions) {
            if (isValidBipartition(bip)) {
                validBipartitions.add(bip);
            }
        }
        
        return validBipartitions;
    }
    
    /**
     * Check if a bipartition is valid
     */
    private boolean isValidBipartition(STBipartition bip) {
        // Check that clusters are non-empty
        if (bip.cluster1.cardinality() == 0 || bip.cluster2.cardinality() == 0) {
            return false;
        }
        
        // Check that clusters don't overlap
        BitSet intersection = (BitSet) bip.cluster1.clone();
        intersection.and(bip.cluster2);
        if (intersection.cardinality() > 0) {
            return false;
        }
        
        // Check that union covers all taxa (or is a proper subset)
        BitSet union = (BitSet) bip.cluster1.clone();
        union.or(bip.cluster2);
        if (union.cardinality() > geneTrees.realTaxaCount) {
            return false;
        }
        
        // Check that it's not a trivial bipartition
        if (bip.cluster1.cardinality() == 1 || bip.cluster2.cardinality() == 1) {
            // Allow trivial bipartitions for now, but could be filtered out
        }
        
        return true;
    }
    
    /**
     * Print expansion summary statistics
     */
    private void printExpansionSummary(int originalCount, int finalCount) {
        System.out.println("\n=== BIPARTITION EXPANSION SUMMARY ===");
        System.out.println("Original bipartitions from gene trees: " + originalCount);
        System.out.println("Distance-based bipartitions added:     " + distanceBasedCount);
        System.out.println("Consensus-based bipartitions added:    " + consensusBasedCount);
        System.out.println("Polytomy resolution bipartitions added: " + polytomyResolutionCount);
        System.out.println("--------------------------------------");
        System.out.println("Total additional bipartitions:         " + (distanceBasedCount + consensusBasedCount + polytomyResolutionCount));
        System.out.println("Final candidate count:                 " + finalCount);
        System.out.println("Expansion factor:                      " + String.format("%.2fx", (double) finalCount / originalCount));
        
        if (finalCount > originalCount) {
            double improvement = ((double) (finalCount - originalCount) / originalCount) * 100;
            System.out.println("Improvement:                           +" + String.format("%.1f%%", improvement) + " more candidates");
        }
    }
    
    /**
     * Get expansion statistics
     */
    public ExpansionStatistics getStatistics() {
        return new ExpansionStatistics(
            originalBipartitionCount,
            distanceBasedCount,
            consensusBasedCount,
            polytomyResolutionCount
        );
    }
    
    /**
     * Statistics class for expansion results
     */
    public static class ExpansionStatistics {
        public final int originalCount;
        public final int distanceBasedCount;
        public final int consensusBasedCount;
        public final int polytomyResolutionCount;
        public final int totalAdded;
        
        public ExpansionStatistics(int originalCount, int distanceBasedCount, 
                                 int consensusBasedCount, int polytomyResolutionCount) {
            this.originalCount = originalCount;
            this.distanceBasedCount = distanceBasedCount;
            this.consensusBasedCount = consensusBasedCount;
            this.polytomyResolutionCount = polytomyResolutionCount;
            this.totalAdded = distanceBasedCount + consensusBasedCount + polytomyResolutionCount;
        }
        
        public double getExpansionFactor() {
            return originalCount > 0 ? (double) (originalCount + totalAdded) / originalCount : 1.0;
        }
    }
}
