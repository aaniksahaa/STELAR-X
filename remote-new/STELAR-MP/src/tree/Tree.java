package tree;

import java.util.ArrayList;
import java.util.Map;
import java.util.Stack;

import taxon.Taxon;

public class Tree {
    
    public ArrayList<TreeNode> nodes;
    public TreeNode root;
    public Map<String, Taxon> taxaMap;
    public TreeNode[] leaves;
    public int leavesCount;
    public boolean isRooted = true; // Always rooted
    
    public TreeNode addNode(ArrayList<TreeNode> children, TreeNode parent) {
        TreeNode nd = new TreeNode().setIndex(nodes.size()).setChilds(children).setParent(parent);
        nodes.add(nd);
        return nd;
    }
    
    public TreeNode addInternalNode(ArrayList<TreeNode> children) {
        var nd = addNode(children, null);
        for (var x : children)
            x.setParent(nd);
        return nd;
    }
    
    public TreeNode addLeaf(Taxon taxon) {
        var nd = addNode(null, null).setTaxon(taxon);
        return nd;
    }
    
    private void parseFromNewick(String newickLine) {
        nodes = new ArrayList<>();
        Stack<TreeNode> nodeStack = new Stack<>();
        newickLine = newickLine.replaceAll("\\s", "");
        
        int n = newickLine.length();
        int i = 0;
        int leavesCount = 0;
        
        while (i < n) {
            char curr = newickLine.charAt(i);
            if (curr == '(') {
                nodeStack.push(null);
                i++;
            } else if (curr == ')') {
                ArrayList<TreeNode> children = new ArrayList<>();
                while (!nodeStack.isEmpty() && nodeStack.peek() != null) {
                    children.add(nodeStack.pop());
                }
                if (!nodeStack.isEmpty())
                    nodeStack.pop();
                
                TreeNode internalNode = addInternalNode(children);
                
                // Skip any support values or branch lengths after ')'
                i++;
                i = skipBranchInfo(newickLine, i);
                
                nodeStack.push(internalNode);
            } else if (curr == ',' || curr == ';') {
                i++;
            } else {
                // Parse taxon name
                StringBuilder taxaName = new StringBuilder();
                while (i < n) {
                    char c = newickLine.charAt(i);
                    if (c == ')' || c == ',' || c == ':' || c == ';') {
                        break;
                    }
                    taxaName.append(c);
                    i++;
                }
                
                // Create leaf node
                String taxonLabel = taxaName.toString();
                if (!taxonLabel.isEmpty()) {
                    Taxon taxon = this.taxaMap.get(taxonLabel);
                    if (taxon == null) {
                        throw new RuntimeException("Unknown taxon: " + taxonLabel);
                    }
                    TreeNode newNode = addLeaf(taxon);
                    leavesCount++;
                    nodeStack.push(newNode);
                }
                
                // Skip any branch length information
                i = skipBranchInfo(newickLine, i);
            }
        }
        
        this.leavesCount = leavesCount;
        if (!nodeStack.isEmpty()) {
            root = nodeStack.lastElement();
        }
        
        // Validate that this is a proper rooted tree
        validateRootedTree();
        
        filterLeaves();
    }
    
    private int skipBranchInfo(String newickLine, int startIndex) {
        int i = startIndex;
        int n = newickLine.length();
        
        // Skip branch length (starts with ':')
        if (i < n && newickLine.charAt(i) == ':') {
            i++; // skip ':'
            // Skip the numeric value (including scientific notation)
            while (i < n) {
                char c = newickLine.charAt(i);
                if (c == ',' || c == ')' || c == ';') {
                    break;
                }
                i++;
            }
        }
        
        return i;
    }
    
    private void validateRootedTree() {
        if (root == null) {
            throw new RuntimeException("Invalid tree: no root found");
        }
        
        if (root.isLeaf()) {
            if (leavesCount != 1) {
                throw new RuntimeException("Invalid tree: single leaf but multiple taxa");
            }
            return; // Single taxon tree is valid
        }
        
        // Root should have exactly 2 children for a proper rooted binary tree
        if (root.childs == null || root.childs.size() != 2) {
            throw new RuntimeException("Invalid rooted tree: root must have exactly 2 children, found: " + 
                (root.childs == null ? 0 : root.childs.size()));
        }
        
        // Validate all internal nodes have exactly 2 children
        validateBinaryStructure(root);
    }
    
    private void validateBinaryStructure(TreeNode node) {
        if (node.isLeaf()) {
            return;
        }
        
        if (node.childs == null || node.childs.size() != 2) {
            throw new RuntimeException("Invalid rooted tree: internal node must have exactly 2 children, found: " + 
                (node.childs == null ? 0 : node.childs.size()));
        }
        
        for (TreeNode child : node.childs) {
            validateBinaryStructure(child);
        }
    }
    
    private void filterLeaves() {
        this.leaves = new TreeNode[this.taxaMap.size()];
        for (var x : nodes) {
            if (x.isLeaf()) {
                this.leaves[x.taxon.id] = x;
            }
        }
    }
    
    public Tree(String newickLine, Map<String, Taxon> taxaMap) {
        this.taxaMap = taxaMap;
        parseFromNewick(newickLine);
    }
    
    public Tree() {
        nodes = new ArrayList<>();
        isRooted = true;
    }
    
    private String newickFormatUtil(TreeNode node) {
        if (node.isLeaf()) {
            return node.taxon.label;
        }
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        for (int i = 0; i < node.childs.size(); ++i) {
            sb.append(newickFormatUtil(node.childs.get(i)));
            if (i != node.childs.size() - 1)
                sb.append(",");
        }
        sb.append(")");
        return sb.toString();
    }
    
    public String getNewickFormat() {
        if (root == null) return "";
        return newickFormatUtil(root) + ";";
    }
    
    public boolean isTaxonPresent(int id) {
        return this.leaves[id] != null;
    }
    
    public boolean isRooted() {
        return true; // Always rooted
    }
    
    // Simplified frequency calculation for rooted trees
    public void calculateFrequencies(Map<String, TreeNode> triPartitions) {
        // Simple implementation - can be enhanced if needed
        calculateFrequenciesUtil(root, triPartitions);
    }
    
    private ArrayList<Integer> calculateFrequenciesUtil(TreeNode node, Map<String, TreeNode> triPartitions) {
        if (node.isLeaf()) {
            ArrayList<Integer> result = new ArrayList<>();
            result.add(node.taxon.id);
            return result;
        }
        
        ArrayList<ArrayList<Integer>> childTaxa = new ArrayList<>();
        for (TreeNode child : node.childs) {
            childTaxa.add(calculateFrequenciesUtil(child, triPartitions));
        }
        
        // Create partition key
        StringBuilder sb = new StringBuilder();
        for (var taxa : childTaxa) {
            taxa.sort(Integer::compareTo);
            for (int taxon : taxa) {
                sb.append(taxon).append("-");
            }
            sb.append("|");
        }
        
        String key = sb.toString();
        if (triPartitions.containsKey(key)) {
            triPartitions.get(key).frequency++;
            node.frequency = 0;
        } else {
            node.frequency = 1;
            triPartitions.put(key, node);
        }
        
        ArrayList<Integer> result = new ArrayList<>();
        for (var taxa : childTaxa) {
            result.addAll(taxa);
        }
        return result;
    }
    
    // Structural equality check (ignores child order)
    public boolean equalsStructure(Tree other) {
        if (other == null || this.root == null || other.root == null) return false;
        return equalsStructureHelper(this.root, other.root);
    }

    private boolean equalsStructureHelper(TreeNode n1, TreeNode n2) {
        if (n1.isLeaf() && n2.isLeaf()) {
            return n1.taxon.label.equals(n2.taxon.label);
        }
        if (n1.isLeaf() || n2.isLeaf()) {
            return false;
        }
        // Both are internal nodes with two children (rooted binary trees)
        if (n1.childs.size() != 2 || n2.childs.size() != 2) return false;
        // Compare both possible child orderings
        TreeNode a1 = n1.childs.get(0), b1 = n1.childs.get(1);
        TreeNode a2 = n2.childs.get(0), b2 = n2.childs.get(1);
        boolean direct = equalsStructureHelper(a1, a2) && equalsStructureHelper(b1, b2);
        boolean swapped = equalsStructureHelper(a1, b2) && equalsStructureHelper(b1, a2);
        return direct || swapped;
    }
}
