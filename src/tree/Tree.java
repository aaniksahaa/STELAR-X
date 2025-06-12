package tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.Stack;

import taxon.Taxon;

public class Tree {
    
    // Core tree structure
    public ArrayList<TreeNode> nodes;               
    public ArrayList<TreeNode> topSortedNodes;      
    
    public TreeNode root;                           // Root node of the tree
    public Map<String, Taxon> taxaMap;          // Mapping from taxon names to RealTaxon objects
    
    // Leaf access optimization for Algorithm 2 scoring
    public TreeNode[] leaves;                       // Fast access to leaf nodes by taxon ID
    public int leavesCount;                         // Number of leaves in tree
    
    public boolean isRooted;

    /**
     * Creates a new internal or leaf tree node.
     * 
     * This is the fundamental building block for tree construction, used both
     * during Newick parsing and during the conquer phase when combining
     * subproblem solutions.
     */
    public TreeNode addNode(ArrayList<TreeNode> children, TreeNode parent){
        TreeNode nd = new TreeNode().setIndex(nodes.size()).setChilds(children).setParent(parent);
        nodes.add(nd);
        return nd;
    }


    /**
     * Creates a new internal node with specified children.
     * 
     * Used extensively during tree construction and the conquer phase where
     * subproblem solutions are combined through tree grafting operations.
     */
    public TreeNode addInternalNode(ArrayList<TreeNode> children){
        var nd = addNode(children, null);
        for (var x : children)
            x.setParent(nd);
        return nd;
    }

    /**
     * Creates a new leaf node for a real taxon.
     * 
     * Leaf nodes represent the actual species/taxa being analyzed. During
     * Algorithm 2 scoring, these nodes are the evaluation points for the
     * mathematical formulations from Section 2.5.
     */
    public TreeNode addLeaf(Taxon taxon){
        var nd = addNode(null, null).setTaxon(taxon);
        return nd;
    }


    /**
     * Parses a phylogenetic tree from Newick format.
     * 
     * This method handles the input gene trees and consensus trees used
     * throughout the wQFM-TREE algorithm. The parsed trees are used for:
     * 1. Gene tree preprocessing and scoring (Algorithm 2)
     * 2. Consensus tree construction for Algorithm 1
     * 3. Final species tree output
     * 
     * The parser handles standard Newick format with taxon name mapping
     * through the provided taxaMap for consistent taxon identification.
     */
    private void parseFromNewick(String newickLine){
        int leavesCount = 0;
        nodes = new ArrayList<>();
        Stack<TreeNode> nodeStack = new Stack<>();
        newickLine = newickLine.replaceAll("\\s", "");
        
        int n = newickLine.length();
        int i = 0, j = 0;
        
        while(i < n){
            char curr = newickLine.charAt(i);
            if(curr == '('){
                nodeStack.push(null);
            }
            else if(curr == ')'){
                ArrayList<TreeNode> children = new ArrayList<>();
                while(!nodeStack.isEmpty() && nodeStack.peek() != null){
                    children.add(nodeStack.pop());
                }
                if(!nodeStack.isEmpty())
                    nodeStack.pop();
                nodeStack.push(addInternalNode(children));
            }
            else if(curr == ',' || curr == ';'){
                // Skip separators
            }
            else{
                StringBuilder taxa = new StringBuilder();
                j = i;
                TreeNode newNode = null;
                while(j < n){
                    char curr_j = newickLine.charAt(j);
                    if(curr_j == ')' || curr_j == ','){
                        Taxon taxon = this.taxaMap.get(taxa.toString());
                        newNode = addLeaf(taxon);
                        leavesCount++;
                        break;
                    }
                    taxa.append(curr_j);
                    ++j;
                }
                if(j == n){
                    leavesCount++;
                    Taxon taxon = this.taxaMap.get(taxa.toString());
                    newNode = addLeaf(taxon);
                }
                i = j - 1;
                nodeStack.push(newNode);
            }
            ++i;
        }

        this.leavesCount = leavesCount;
        root = nodeStack.lastElement();
        
        // Determine if tree is rooted or unrooted
        determineRootedness();
        
        // Handle tree structure based on rootedness
        if(!isRooted && root.childs != null && root.childs.size() > 2){
            balanceRoot();
        }
        
        filterLeaves();
        topSort();
    }

    private void determineRootedness(){
        if(root.childs == null){
            isRooted = true; // Single leaf case
            return;
        }
        
        // Unrooted trees typically have 3 children at root
        // Rooted trees typically have 2 children at root
        isRooted = (root.childs.size() <= 2);
    }

    public void setRooted(boolean rooted){
        this.isRooted = rooted;
        if(!isRooted && root.childs != null && root.childs.size() > 2){
            balanceRoot();
        }
    }

    /**
     * Creates fast-access array for leaf nodes indexed by taxon ID.
     * 
     * This optimization is crucial for Algorithm 2 scoring, which needs
     * frequent access to leaf nodes by taxon ID during quartet evaluation.
     * The leaves array enables O(1) lookup instead of O(n) tree traversal.
     */
    private void filterLeaves(){
        this.leaves = new TreeNode[this.taxaMap.size()];
        for(var x : nodes){
            if(x.isLeaf()){
                this.leaves[x.taxon.id] = x;
            }
        }
    }

    /**
     * Recursively resolves non-binary internal nodes using distance matrix.
     * 
     * This method converts non-binary (polytomy) nodes into binary nodes using
     * a distance-based heuristic. This is important for Algorithm 2 scoring
     * which assumes binary tree structure for efficient quartet evaluation.
     * 
     * The resolution strategy selects the child subtree with maximum total
     * distance to all other taxa and separates it from the rest.
     */
    public ArrayList<Integer> resolveNonBinaryUtil(TreeNode node, double[][] distanceMatrix){
        if(node.isLeaf()){
            return new ArrayList<>(Arrays.asList(node.taxon.id));
        }
        
        if(node.childs.size() <= 2){
            ArrayList<Integer> result = new ArrayList<>();
            for(var child : node.childs){
                result.addAll(resolveNonBinaryUtil(child, distanceMatrix));
            }
            return result;
        }

        double maxDistance = -1;
        int maxIndex = -1;
        
        for(int i = 0; i < node.childs.size(); i++){
            var childTaxa = resolveNonBinaryUtil(node.childs.get(i), distanceMatrix);
            double totalDistance = 0;
            
            for(int taxon1 : childTaxa){
                for(int j = 0; j < node.childs.size(); j++){
                    if(i == j) continue;
                    var otherChildTaxa = resolveNonBinaryUtil(node.childs.get(j), distanceMatrix);
                    for(int taxon2 : otherChildTaxa){
                        totalDistance += distanceMatrix[taxon1][taxon2];
                    }
                }
            }
            
            if(totalDistance > maxDistance){
                maxDistance = totalDistance;
                maxIndex = i;
            }
        }

        TreeNode maxChild = node.childs.get(maxIndex);
        ArrayList<TreeNode> remainingChildren = new ArrayList<>();
        for(int i = 0; i < node.childs.size(); i++){
            if(i != maxIndex){
                remainingChildren.add(node.childs.get(i));
            }
        }

        TreeNode newInternalNode = addInternalNode(remainingChildren);
        
        ArrayList<TreeNode> newChildren = new ArrayList<>();
        newChildren.add(maxChild);
        newChildren.add(newInternalNode);
        
        node.childs = newChildren;
        maxChild.parent = node;
        newInternalNode.parent = node;

        ArrayList<Integer> result = new ArrayList<>();
        result.addAll(resolveNonBinaryUtil(maxChild, distanceMatrix));
        result.addAll(resolveNonBinaryUtil(newInternalNode, distanceMatrix));
        return result;
    }

    public void resolveNonBinary(double[][] distanceMatrix){
        for(var node : nodes){
            if(node.childs != null && node.childs.size() > 2){
                resolveNonBinaryUtil(node, distanceMatrix);
            }
        }
    }

    public Tree(String newickLine, Map<String, Taxon> taxaMap){
        this.taxaMap = taxaMap;
        parseFromNewick(newickLine);
    }

    public Tree(){
        nodes = new ArrayList<>();
    }

    public int dfs(TreeNode node, ArrayList<Integer> subTreeNodeCount){
        if(node.isLeaf()){
            subTreeNodeCount.set(node.index, 1);
            return 1;
        }
        
        int count = 1;
        for(var child : node.childs){
            count += dfs(child, subTreeNodeCount);
        }
        subTreeNodeCount.set(node.index, count);
        return count;
    }

    public void reRootTree(TreeNode newRootNode){
        TreeNode newRootP = newRootNode.parent;
        if(newRootP == null) return;
        newRootP.childs.remove(newRootNode);

        TreeNode curr = newRootP;
        TreeNode currP, temp;
        currP = curr.parent;
        while(curr != null && currP != null){
            curr.childs.add(currP);
            currP.childs.remove(curr);
            temp = currP;
            currP = currP.parent;
            temp.parent = curr;
            curr = temp;
        }
        if(newRootNode.isLeaf())
            newRootNode.childs = new ArrayList<>();
        newRootNode.childs.add(newRootP);
        this.root = newRootNode;
        this.isRooted = true;
    }

    private void balanceRoot(){
        int n = nodes.size();
        ArrayList<Integer> subTreeNodeCount = new ArrayList<>(n);
        for(int i = 0; i < n; ++i)
            subTreeNodeCount.add(0);
        dfs(root, subTreeNodeCount);
        
        TreeNode closest = root;
        int diff = n;
        int v;
        for(int i = 0; i < n; ++i){
            v = Math.abs(n/2 - subTreeNodeCount.get(i)); 
            if(v < diff){
                diff = v;
                closest = nodes.get(i);
            }
        }
        
        TreeNode closestP = closest.parent;
        if(closestP == null) return;
        
        closestP.childs.remove(closest);

        TreeNode curr = closestP;
        TreeNode currP, temp;
        currP = curr.parent;
        while(curr != null && currP != null){
            curr.childs.add(currP);
            currP.childs.remove(curr);
            temp = currP;
            currP = currP.parent;
            temp.parent = curr;
            curr = temp;
        }

        ArrayList<TreeNode> arr = new ArrayList<>();
        arr.add(closest);
        arr.add(closestP);

        root = addInternalNode(arr);
    }
    
    private String newickFormatUtil(TreeNode node){
        if(node.isLeaf()){
            return node.taxon.label;
        }
        StringBuilder sb = new StringBuilder();
        sb.append("(");
        for(int i = 0; i < node.childs.size(); ++i){
            sb.append(newickFormatUtil(node.childs.get(i)));
            if(i != node.childs.size() - 1)
                sb.append(",");
        }
        sb.append(")");
        return sb.toString();
    }

    public String getNewickFormat(){
        return newickFormatUtil(root) + ";";
    }

    public boolean isTaxonPresent(int id){
        return this.leaves[id] != null;
    }
    
    private ArrayList<Integer> getChildrens(TreeNode node, Map<String, TreeNode> triPartitionsMap){
        if(node.isLeaf()){
            ArrayList<Integer> arr = new ArrayList<>();
            arr.add(node.taxon.id);
            return arr;
        }
        ArrayList<ArrayList<Integer>> reachableFromChilds = new ArrayList<>();
        for(var child : node.childs){
            reachableFromChilds.add(getChildrens(child, triPartitionsMap));
        }

        boolean[] mark = new boolean[this.taxaMap.size()];
        for(var childTaxa : reachableFromChilds){
            for(var x : childTaxa){
                mark[x] = true;
            }
        }
        
        ArrayList<Integer> arr = new ArrayList<>();
        for(int i = 0; i < this.taxaMap.size(); ++i){
            if(!mark[i] && isTaxonPresent(i)){
                arr.add(i);
            }
        }

        reachableFromChilds.add(arr);
        String[] partitionStrings = new String[reachableFromChilds.size()];

        for(var x : reachableFromChilds){
            x.sort((a, b) -> a - b);
        }

        for(int i = 0; i < reachableFromChilds.size(); ++i){
            var sb = new StringBuilder();
            reachableFromChilds.get(i).forEach(s -> sb.append(s + '-') );
            partitionStrings[i] = sb.toString();
        }
        
        Arrays.sort(partitionStrings);
        StringBuilder sb = new StringBuilder();
        for(var x : partitionStrings){
            sb.append(x + '|');
        }

        String key = sb.toString();

        if(triPartitionsMap.containsKey(key)){
            triPartitionsMap.get(key).frequency++;
            node.frequency = 0;
        }
        else{
            node.frequency = 1;
            triPartitionsMap.put(key, node);
        }
        reachableFromChilds.remove(reachableFromChilds.size() - 1);
        
        arr = new ArrayList<>();
        for(var childTaxa : reachableFromChilds){
            arr.addAll(childTaxa);
        }
        return arr;
    }

    public void calculateFrequencies(Map<String, TreeNode> triPartitions){
        getChildrens(root, triPartitions);
    }

    private void topSortUtil(TreeNode node, ArrayList<TreeNode> topSort){
        if(!node.isLeaf()){
            for(var x : node.childs){
                topSortUtil(x, topSort);
            }
        }
        topSort.add(node);
    }

    public void topSort(){
        ArrayList<TreeNode> topSort = new ArrayList<>();
        topSortUtil(root, topSort);
        this.topSortedNodes = topSort;
    }

    public boolean checkIfNonBinary(){
        for(var x : nodes){
            if(x.childs != null && x.childs.size() > 2)
                return true;
        }
        return false;
    }

    public boolean isRooted(){
        return isRooted;
    }

    public void forceUnrooted(){
        this.isRooted = false;
        if(root.childs != null && root.childs.size() == 2){
            // Convert rooted to unrooted by adding artificial root
            TreeNode oldRoot = root;
            ArrayList<TreeNode> newRootChildren = new ArrayList<>();
            newRootChildren.addAll(oldRoot.childs);
            root = addInternalNode(newRootChildren);
        }
    }

    public void forceRooted(){
        this.isRooted = true;
        if(root.childs != null && root.childs.size() > 2){
            balanceRoot();
        }
    }
}
