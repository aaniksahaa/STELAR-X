// Debug script to compare original vs new bipartition generation
// Original algorithm that we need to replicate:

// For each internal node in a rooted binary tree:
// 1. Calculate BitSets for left and right children recursively
// 2. Create STBipartition(leftCluster, rightCluster) 
// 3. Return union of left and right clusters for this node

// Key insight: We get one bipartition per internal node
// For n taxa, we get n-1 internal nodes, so n-1 bipartitions per tree

// The issue might be:
// 1. We're only extracting left subtrees (missing right subtrees)
// 2. Range calculation is incorrect  
// 3. We're not handling the tree structure properly
