    }
    
    // Removed expensive buildCandidateRangeMapping and findMatchingRange methods
    // Using direct calculation approach instead
    
    /**
     * Convert BitSet to set of taxon IDs.
     */
    private Set<Integer> bitSetToTaxonSet(utils.BitSet bitSet) {
        Set<Integer> taxonSet = new HashSet<>();
        for (int i = bitSet.nextSetBit(0); i >= 0; i = bitSet.nextSetBit(i + 1)) {
            taxonSet.add(i);
        }
        return taxonSet;
    }
    
    /**
     * Get taxon set for a range in a gene tree.
     */
    private Set<Integer> getRangeTaxonSet(int geneTreeIndex, int start, int end) {
        Set<Integer> taxonSet = new HashSet<>();
        int[][] orderings = bipartitionManager.getGeneTreeTaxaOrdering();
        
        if (orderings != null && geneTreeIndex < orderings.length) {
            int[] ordering = orderings[geneTreeIndex];
            if (ordering != null) {
                for (int pos = start; pos < end && pos < ordering.length; pos++) {
                    if (ordering[pos] != -1) { // Valid taxon
                        taxonSet.add(ordering[pos]);
                    }
                }
            }
        }
        
        return taxonSet;
    }
}
