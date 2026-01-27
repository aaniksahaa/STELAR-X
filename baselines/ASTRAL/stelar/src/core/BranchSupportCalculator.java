package core;

import java.util.Map;
import java.text.DecimalFormat;

import preprocessing.GeneTrees;
import tree.Tree;
import tree.TreeNode;
import core.QuartetFrequencyCounter.QuartetFrequencies;

/**
 * Main class for calculating branch support and branch lengths for species trees.
 * 
 * This class integrates quartet frequency counting with posterior probability calculation
 * to provide ASTRAL-style branch annotations including:
 * - Posterior probabilities
 * - Branch lengths based on coalescent theory
 * - P-values for polytomy testing
 * 
 * Usage:
 * 1. Create instance with gene trees and inferred species tree
 * 2. Call annotateBranches() to add support values
 * 3. Access annotated tree with support values
 */
public class BranchSupportCalculator {
    
    private GeneTrees geneTrees;
    private Tree speciesTree;
    private double lambda; // Prior parameter for Beta-Binomial model
    private BranchAnnotationType annotationType;
    
    private DecimalFormat df = new DecimalFormat("#.####");
    
    /**
     * Types of branch annotations supported
     */
    public enum BranchAnnotationType {
        NONE,                    // No annotation
        POSTERIOR_ONLY,          // Only posterior probabilities
        DETAILED,               // Posterior, frequencies, effective N
        BRANCH_LENGTH_ONLY,     // Only branch lengths
        POSTERIOR_AND_LENGTH,   // Both posterior and branch length
        PVALUE_ONLY,            // Only p-values for polytomy test
        ALL                     // All available information
    }
    
    /**
     * Constructor
     * 
     * @param geneTrees The input gene trees
     * @param speciesTree The inferred species tree to annotate
     * @param lambda Prior parameter (typically 0.5)
     * @param annotationType Type of annotation to perform
     */
    public BranchSupportCalculator(GeneTrees geneTrees, Tree speciesTree, 
                                 double lambda, BranchAnnotationType annotationType) {
        this.geneTrees = geneTrees;
        this.speciesTree = speciesTree;
        this.lambda = lambda;
        this.annotationType = annotationType;
    }
    
    /**
     * Constructor with default parameters
     */
    public BranchSupportCalculator(GeneTrees geneTrees, Tree speciesTree) {
        this(geneTrees, speciesTree, 0.5, BranchAnnotationType.POSTERIOR_ONLY);
    }
    
    /**
     * Main method to annotate all branches in the species tree with support values
     */
    public void annotateBranches() {
        if (annotationType == BranchAnnotationType.NONE) {
            return;
        }
        
        System.out.println("Calculating branch support using quartet frequencies...");
        
        // Count quartet frequencies for all branches
        QuartetFrequencyCounter counter = new QuartetFrequencyCounter(geneTrees, speciesTree);
        Map<TreeNode, QuartetFrequencies> branchFrequencies = counter.countAllBranchFrequencies();
        
        System.out.println("Found " + branchFrequencies.size() + " internal branches to annotate");
        
        // Calculate posterior probabilities and branch lengths
        int annotatedBranches = 0;
        double totalLogPosterior = 0.0;
        
        for (Map.Entry<TreeNode, QuartetFrequencies> entry : branchFrequencies.entrySet()) {
            TreeNode node = entry.getKey();
            QuartetFrequencies freqs = entry.getValue();
            
            // Skip branches with insufficient data
            if (freqs.n < 1) {
                node.supportValue = null;
                continue;
            }
            
            // Calculate posterior probabilities and branch length
            PosteriorCalculator posterior = new PosteriorCalculator(
                freqs.f1, freqs.f2, freqs.f3, freqs.n, lambda);
            
            double postProb = posterior.getPosteriorProbability();
            double branchLength = posterior.getBranchLength();
            
            
            // Set branch length in the tree structure
            setBranchLength(node, branchLength);
            
            // Create annotation string based on type
            String annotation = createAnnotation(posterior, freqs);
            node.supportValue = annotation;
            
            // Accumulate log posterior for total score
            if (postProb > 0) {
                totalLogPosterior += Math.log(postProb);
            }
            
            annotatedBranches++;
            
            // Warn about low effective N
            if (freqs.n < 20) {
                System.err.println("Warning: Branch with low effective N (" + freqs.n + 
                                 ") may have unreliable support values");
            }
            
        }
        
        System.out.println("Annotated " + annotatedBranches + " branches with support values");
        System.out.println("Total log posterior probability: " + df.format(totalLogPosterior));
        
        // Set precision for decimal formatting
        df.setMaximumFractionDigits(4);
    }
    
    /**
     * Set branch length for a node
     */
    private void setBranchLength(TreeNode node, double branchLength) {
        // Store branch length in the node
        // This could be added to TreeNode class if needed
        node.branchLength = branchLength;
    }
    
    /**
     * Create annotation string based on annotation type
     */
    private String createAnnotation(PosteriorCalculator posterior, QuartetFrequencies freqs) {
        switch (annotationType) {
            case POSTERIOR_ONLY:
                return df.format(posterior.getPosteriorProbability());
                
            case BRANCH_LENGTH_ONLY:
                return df.format(posterior.getBranchLength());
                
            case POSTERIOR_AND_LENGTH:
                return df.format(posterior.getPosteriorProbability()) + 
                       ":" + df.format(posterior.getBranchLength());
                
            case PVALUE_ONLY:
                double pval = posterior.getPValue();
                if (pval < 0) {
                    return "NA";
                } else {
                    return df.format(pval);
                }
                
            case DETAILED:
                return String.format("'[pp=%.4f;f1=%.1f;f2=%.1f;f3=%.1f;EN=%.1f]'",
                    posterior.getPosteriorProbability(),
                    freqs.f1, freqs.f2, freqs.f3, freqs.n);
                
            case ALL:
                return String.format("'[pp=%.4f;bl=%.4f;f1=%.1f;f2=%.1f;f3=%.1f;EN=%.1f;pval=%.4f]'",
                    posterior.getPosteriorProbability(),
                    posterior.getBranchLength(),
                    freqs.f1, freqs.f2, freqs.f3, freqs.n,
                    posterior.getPValue() < 0 ? Double.NaN : posterior.getPValue());
                
            default:
                return null;
        }
    }
    
    /**
     * Get the annotated species tree
     */
    public Tree getAnnotatedTree() {
        return speciesTree;
    }
    
    /**
     * Calculate and return summary statistics about branch support
     */
    public BranchSupportStatistics calculateStatistics() {
        QuartetFrequencyCounter counter = new QuartetFrequencyCounter(geneTrees, speciesTree);
        Map<TreeNode, QuartetFrequencies> branchFrequencies = counter.countAllBranchFrequencies();
        
        int totalBranches = branchFrequencies.size();
        int highSupportBranches = 0;  // posterior > 0.95
        int mediumSupportBranches = 0; // posterior 0.75-0.95
        int lowSupportBranches = 0;    // posterior < 0.75
        
        double avgPosterior = 0.0;
        double avgEffectiveN = 0.0;
        double totalLogPosterior = 0.0;
        
        int validBranches = 0;
        
        for (QuartetFrequencies freqs : branchFrequencies.values()) {
            if (freqs.n < 1) continue;
            
            PosteriorCalculator posterior = new PosteriorCalculator(
                freqs.f1, freqs.f2, freqs.f3, freqs.n, lambda);
            
            double postProb = posterior.getPosteriorProbability();
            
            if (postProb >= 0.95) {
                highSupportBranches++;
            } else if (postProb >= 0.75) {
                mediumSupportBranches++;
            } else {
                lowSupportBranches++;
            }
            
            avgPosterior += postProb;
            avgEffectiveN += freqs.n;
            if (postProb > 0) {
                totalLogPosterior += Math.log(postProb);
            }
            
            validBranches++;
        }
        
        if (validBranches > 0) {
            avgPosterior /= validBranches;
            avgEffectiveN /= validBranches;
        }
        
        return new BranchSupportStatistics(
            totalBranches, validBranches,
            highSupportBranches, mediumSupportBranches, lowSupportBranches,
            avgPosterior, avgEffectiveN, totalLogPosterior
        );
    }
    
    /**
     * Statistics about branch support calculation
     */
    public static class BranchSupportStatistics {
        public final int totalBranches;
        public final int validBranches;
        public final int highSupportBranches;  // >= 0.95
        public final int mediumSupportBranches; // 0.75-0.95
        public final int lowSupportBranches;   // < 0.75
        public final double averagePosterior;
        public final double averageEffectiveN;
        public final double totalLogPosterior;
        
        public BranchSupportStatistics(int total, int valid, int high, int medium, int low,
                                     double avgPost, double avgN, double logPost) {
            this.totalBranches = total;
            this.validBranches = valid;
            this.highSupportBranches = high;
            this.mediumSupportBranches = medium;
            this.lowSupportBranches = low;
            this.averagePosterior = avgPost;
            this.averageEffectiveN = avgN;
            this.totalLogPosterior = logPost;
        }
        
        @Override
        public String toString() {
            DecimalFormat df = new DecimalFormat("#.####");
            StringBuilder sb = new StringBuilder();
            sb.append("Branch Support Statistics:\n");
            sb.append("  Total branches: ").append(totalBranches).append("\n");
            sb.append("  Valid branches: ").append(validBranches).append("\n");
            sb.append("  High support (>=0.95): ").append(highSupportBranches).append("\n");
            sb.append("  Medium support (0.75-0.95): ").append(mediumSupportBranches).append("\n");
            sb.append("  Low support (<0.75): ").append(lowSupportBranches).append("\n");
            sb.append("  Average posterior: ").append(df.format(averagePosterior)).append("\n");
            sb.append("  Average effective N: ").append(df.format(averageEffectiveN)).append("\n");
            sb.append("  Total log posterior: ").append(df.format(totalLogPosterior));
            return sb.toString();
        }
    }
    
    /**
     * Validate quartet frequencies for debugging
     */
    public void validateQuartetFrequencies() {
        System.out.println("Validating quartet frequency calculations...");
        
        QuartetFrequencyCounter counter = new QuartetFrequencyCounter(geneTrees, speciesTree);
        Map<TreeNode, QuartetFrequencies> branchFrequencies = counter.countAllBranchFrequencies();
        
        for (Map.Entry<TreeNode, QuartetFrequencies> entry : branchFrequencies.entrySet()) {
            QuartetFrequencies freqs = entry.getValue();
            
            // Check for reasonable frequency values
            if (freqs.f1 < 0 || freqs.f2 < 0 || freqs.f3 < 0) {
                System.err.println("Warning: Negative frequencies detected: " + freqs);
            }
            
            if (Math.abs(freqs.f1 + freqs.f2 + freqs.f3 - freqs.n) > 0.001) {
                System.err.println("Warning: Frequency sum mismatch: " + freqs);
            }
            
            if (freqs.n > freqs.quartCount) {
                System.err.println("Warning: Effective N exceeds total quartets: " + freqs);
            }
        }
        
        System.out.println("Quartet frequency validation completed");
    }
    
}
