package stelarx.tree;

import stelarx.Logging;
import stelarx.taxon.TaxonRegistry;
import stelarx.util.ProgressBar;
import stelarx.util.Threading;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;

/**
 * Newick parser for explicitly rooted gene trees with optional native polytomies.
 *
 * Two-pass design:
 *   Pass 1  -- collectTaxonNames(): scan every Newick string, register all names.
 *   Pass 2  -- parseNewick(): build Tree objects with postorder arrays + node ranges.
 *
 * Parsing is intentionally lenient: any number of children is allowed during the
 * stack-based parse phase.  The conversion step then applies the selected input
 * policy and roots the in-memory representation:
 *
 *   Root node with 2 children → rooted binary tree, keep its supplied root.
 *   Any other root arity → reject: STELAR-X never invents an arbitrary root.
 *   Default                            → deterministic first-pair binary refinement.
 *   keepPolytomy=true internal degree≥3 → native polytomous node.
 *
 * After parsing every node has a half-open range [rangeStart, rangeEnd) that indexes
 * into the tree's postorderArray (left-to-right leaf ordering).
 */
public class TreeParser {

    // -------------------------------------------------------------------------
    // Public entry point
    // -------------------------------------------------------------------------

    public static List<Tree> parseGeneTrees(String inputFile,
                                             TaxonRegistry registry) throws IOException {
        return parseGeneTrees(inputFile, registry, false);
    }

    /**
     * Parse gene trees and, by default, apply the same deterministic first-pair
     * binary refinement available in {@code clean.py --deterministic}. Passing
     * {@code keepPolytomy=true} retains native multifurcations for downstream
     * unresolved rooted-triplet scoring.
     */
    public static List<Tree> parseGeneTrees(String inputFile,
                                             TaxonRegistry registry,
                                             boolean keepPolytomy) throws IOException {
        return parseGeneTreesDetailed(inputFile, registry, keepPolytomy).trees();
    }

    /**
     * Parse gene trees and retain whether genuine unresolved internal
     * multifurcations were observed before any deterministic refinement.
     */
    public static ParsedGeneTrees parseGeneTreesDetailed(String inputFile,
                                                          TaxonRegistry registry,
                                                          boolean keepPolytomy)
            throws IOException {
        long t0 = System.nanoTime();

        // Read all non-empty lines
        List<String> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {
            String ln;
            while ((ln = br.readLine()) != null) {
                ln = ln.trim();
                if (!ln.isEmpty()) lines.add(ln);
            }
        }
        if (lines.isEmpty()) {
            throw new IllegalArgumentException("Tree file is empty: " + inputFile);
        }
        Logging.info("Read %d lines from %s", lines.size(), inputFile);

        // Pass 1 – register taxon names
        for (String ln : lines) collectTaxonNames(ln, registry);
        registry.lock();
        int n = registry.size();
        Logging.info("Registered %d unique taxa", n);

        // Pass 2 is independent after the registry is locked.  In normal CLI runs
        // Threading is already initialized, so parsing + optional refinement uses
        // the configured worker pool while preserving input order in parsed[].
        Tree[] parsed = new Tree[lines.size()];
        int[][] perTreeCounts = new int[lines.size()][6];
        ProgressBar parseBar = new ProgressBar("Parsing trees", lines.size());
        AtomicInteger completed = new AtomicInteger(0);
        if (lines.size() > 1 && Threading.isStarted() && Threading.getNumThreads() > 1) {
            Threading.processRangeParallel(lines.size(), i -> {
                parsed[i] = parseNewick(lines.get(i), i, registry, perTreeCounts[i],
                    keepPolytomy);
                parseBar.update(completed.incrementAndGet());
            });
        } else {
            for (int i = 0; i < lines.size(); i++) {
                parsed[i] = parseNewick(lines.get(i), i, registry, perTreeCounts[i],
                    keepPolytomy);
                parseBar.update(i + 1);
            }
        }
        parseBar.done();
        List<Tree> trees = new ArrayList<>(Arrays.asList(parsed));

        int[] totals = new int[6];
        int treesWithPolytomies = 0;
        int treesRefined = 0;
        for (int[] counts : perTreeCounts) {
            for (int j = 0; j < totals.length; j++) totals[j] += counts[j];
            if (counts[4] > 0) treesWithPolytomies++;
            if (counts[3] > 0) treesRefined++;
        }

        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Parsed %d gene trees in %d ms", trees.size(), ms);
        if (keepPolytomy) {
            Logging.info("Input polytomies: keeping %d node(s) across %d tree(s); "
                + "native downstream handling enabled", totals[4], treesWithPolytomies);
        } else if (totals[3] > 0) {
            Logging.info("Input binary refinement: detected %d polytomy node(s); resolved "
                + "%d multifurcation(s) across %d tree(s) by deterministic first-pair pairing "
                + "(%s)",
                totals[4], totals[3], treesRefined,
                lines.size() > 1 && Threading.isStarted() && Threading.getNumThreads() > 1
                    ? "parallel" : "serial");
        } else {
            Logging.info("Input binary refinement: no multifurcations detected");
        }

        // Per-tree debug log -- cap at 5 trees to avoid flooding on large inputs
        if (Logging.isDebug()) {
            int cap = Math.min(5, trees.size());
            for (int i = 0; i < cap; i++) {
                Tree t = trees.get(i);
                Logging.debug("  Tree %d: %d leaves, complete=%b  postorder=%s",
                    i, t.leafCount, t.isComplete,
                    Logging.isTrace() ? Arrays.toString(t.postorderArray) : "(use -vvv to see)");
            }
            if (trees.size() > cap)
                Logging.debug("  ... (%d more trees not shown)", trees.size() - cap);
        }

        return new ParsedGeneTrees(trees, totals[4]);
    }

    public record ParsedGeneTrees(List<Tree> trees,
                                  int detectedPolytomyNodeCount) {}

    /**
     * Parse induced gene trees directly against a precomputed, locked taxon
     * registry. Leaves outside the registry are discarded while the temporary
     * Newick topology is assembled, before TreeNode objects or compact arrays
     * are created. Unary nodes produced by filtering are suppressed.
     *
     * Trees retaining zero or one selected taxon are discarded because they
     * cannot contribute a rooted triplet. Two-taxon trees are retained at this
     * parsing layer because their clades may still be useful to callers; the
     * inference restriction layer drops them before scoring.
     */
    public static RestrictedGeneTrees parseRestrictedGeneTrees(
            String inputFile, TaxonRegistry registry, boolean keepPolytomy)
            throws IOException {
        if (!registry.isLocked()) {
            throw new IllegalArgumentException(
                "Restricted gene-tree parsing requires a locked taxon registry");
        }

        long t0 = System.nanoTime();
        List<String> lines = readNonEmptyLines(inputFile);
        if (lines.isEmpty()) {
            throw new IllegalArgumentException("Tree file is empty: " + inputFile);
        }
        Logging.info("Read %d lines from %s", lines.size(), inputFile);
        Logging.info("Registered %d selected taxa", registry.size());

        Tree[] parsed = new Tree[lines.size()];
        int[][] perTreeCounts = new int[lines.size()][6];
        ProgressBar parseBar = new ProgressBar("Parsing filtered trees", lines.size());
        AtomicInteger completed = new AtomicInteger(0);
        if (lines.size() > 1 && Threading.isStarted() && Threading.getNumThreads() > 1) {
            Threading.processRangeParallel(lines.size(), i -> {
                parsed[i] = parseNewickRestricted(lines.get(i), i, registry,
                    perTreeCounts[i], keepPolytomy);
                parseBar.update(completed.incrementAndGet());
            });
        } else {
            for (int i = 0; i < lines.size(); i++) {
                parsed[i] = parseNewickRestricted(lines.get(i), i, registry,
                    perTreeCounts[i], keepPolytomy);
                parseBar.update(i + 1);
            }
        }
        parseBar.done();

        List<Tree> trees = new ArrayList<>(lines.size());
        int dropped = 0;
        for (Tree tree : parsed) {
            if (tree == null) {
                dropped++;
                continue;
            }
            int treeIndex = trees.size();
            trees.add(tree.treeIndex == treeIndex ? tree
                : new Tree(treeIndex, tree.root, tree.postorderArray, tree.positionMap,
                    tree.leafCount, registry.size(), tree.hasPolytomy));
        }
        if (trees.isEmpty()) {
            throw new IllegalArgumentException(
                "No gene tree retains at least two taxa after applying the taxa file");
        }

        int[] totals = new int[6];
        for (int[] counts : perTreeCounts) {
            for (int j = 0; j < totals.length; j++) totals[j] += counts[j];
        }
        long ms = (System.nanoTime() - t0) / 1_000_000;
        Logging.info("Parsed %d induced gene trees in %d ms; discarded %d tree(s) "
                + "with fewer than two selected taxa",
            trees.size(), ms, dropped);
        if (keepPolytomy) {
            Logging.info("Filtered input polytomies: keeping %d node(s)", totals[4]);
        } else if (totals[3] > 0) {
            Logging.info("Filtered input binary refinement: resolved %d induced "
                    + "multifurcation(s) after taxon restriction",
                totals[3]);
        }

        if (Logging.isDebug()) {
            int cap = Math.min(5, trees.size());
            for (int i = 0; i < cap; i++) {
                Tree tree = trees.get(i);
                Logging.debug("  Induced tree %d: %d leaves, complete=%b  postorder=%s",
                    i, tree.leafCount, tree.isComplete,
                    Logging.isTrace() ? Arrays.toString(tree.postorderArray)
                        : "(use -vvv to see)");
            }
            if (trees.size() > cap) {
                Logging.debug("  ... (%d more trees not shown)", trees.size() - cap);
            }
        }
        return new RestrictedGeneTrees(trees, lines.size(), dropped, totals[4]);
    }

    public record RestrictedGeneTrees(List<Tree> trees,
                                      int sourceTreeCount,
                                      int droppedTreeCount,
                                      int detectedPolytomyNodeCount) {}

    /**
     * Parse one supplied species tree against an already-locked gene-tree taxon
     * registry. The tree must contain exactly the same taxa as the gene-tree
     * input; unknown, missing, or duplicate taxa are rejected.
     */
    public static Tree parseSpeciesTree(String inputFile,
                                        TaxonRegistry registry) throws IOException {
        if (!registry.isLocked()) {
            throw new IllegalArgumentException("Species-tree parsing requires a locked taxon registry");
        }

        List<String> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {
            String ln;
            while ((ln = br.readLine()) != null) {
                ln = ln.trim();
                if (!ln.isEmpty()) lines.add(ln);
            }
        }
        if (lines.isEmpty()) {
            throw new IllegalArgumentException("Species tree file is empty: " + inputFile);
        }
        if (lines.size() != 1) {
            throw new IllegalArgumentException("Species tree file must contain exactly one Newick tree: " + inputFile);
        }

        int[] rootingCounts = new int[6];
        Tree tree = parseNewick(lines.get(0), 0, registry, rootingCounts, true);
        validateCompleteTaxonSet(tree, registry, inputFile);
        Logging.info("Parsed supplied species tree: %d leaves", tree.leafCount);
        return tree;
    }

    /** Parse an in-memory inferred species-tree Newick against a locked registry. */
    public static Tree parseSpeciesTreeNewick(String newick,
                                              TaxonRegistry registry) {
        if (!registry.isLocked()) {
            throw new IllegalArgumentException(
                "Species-tree parsing requires a locked taxon registry");
        }
        if (newick == null || newick.isBlank()) {
            throw new IllegalArgumentException("Species-tree Newick is empty");
        }

        int[] rootingCounts = new int[6];
        Tree tree = parseNewick(newick.trim(), 0, registry, rootingCounts, true);
        validateCompleteTaxonSet(tree, registry, "inferred species tree");
        return tree;
    }

    /**
     * Parse exactly one tree with its own registry, without requiring it to have
     * the same taxon set as a previously parsed gene-tree collection.  This is
     * used only by opt-in taxon-restricted scoring; the ordinary strict species-
     * tree parser above remains unchanged.
     */
    public static StandaloneTree parseStandaloneTree(String inputFile) throws IOException {
        List<String> lines = readNonEmptyLines(inputFile);
        if (lines.isEmpty()) {
            throw new IllegalArgumentException("Tree file is empty: " + inputFile);
        }
        if (lines.size() != 1) {
            throw new IllegalArgumentException(
                "Expected exactly one Newick tree: " + inputFile);
        }

        TaxonRegistry registry = new TaxonRegistry();
        collectTaxonNames(lines.get(0), registry);
        registry.lock();
        int[] rootingCounts = new int[6];
        Tree tree = parseNewick(lines.get(0), 0, registry, rootingCounts, true);
        return new StandaloneTree(tree, registry);
    }

    public record StandaloneTree(Tree tree, TaxonRegistry registry) {}

    // -------------------------------------------------------------------------
    // Pass 1 – name collection
    // -------------------------------------------------------------------------

    /**
     * Walk the Newick string and register every taxon name.
     * Tracks whether the last structural token was ')': if so, the next label is
     * a bootstrap/internal value (skip it); otherwise it is a taxon name (register it).
     * This correctly handles both named taxa (strings) and integer-labelled taxa.
     */
    private static void collectTaxonNames(String s, TaxonRegistry reg) {
        forEachTaxonName(s, reg::register);
    }

    /**
     * Shared Newick leaf-token scanner.  Taxa extraction and coverage reporting
     * deliberately use this exact scanner so their name semantics cannot drift
     * from normal STELAR-X parsing.
     */
    static void forEachTaxonName(String s, Consumer<String> consumer) {
        int i = 0, n = s.length();
        // true if the most recent structural character was ')'
        boolean afterCloseParen = false;
        while (i < n) {
            char c = s.charAt(i);
            if (c == '(') { afterCloseParen = false; i++; continue; }
            if (c == ',') { afterCloseParen = false; i++; continue; }
            if (c == ')') { afterCloseParen = true;  i++; continue; }
            if (c == ';') break;
            if (c == ':') { i = skipBranchLen(s, i + 1, n); continue; }
            if (c == '[') { // NHX or comment: skip to ']'
                while (i < n && s.charAt(i) != ']') i++;
                if (i < n) i++;
                continue;
            }
            // Token: taxon name (after '(' or ',') or internal label (after ')')
            int start = i;
            i = scanLabelEnd(s, i, n);
            String tok = decodeNewickLabel(s.substring(start, i));
            if (!tok.isEmpty() && !afterCloseParen) {
                consumer.accept(tok);
            }
            afterCloseParen = false;
        }
    }

    static List<String> readNonEmptyLines(String inputFile) throws IOException {
        List<String> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) lines.add(line);
            }
        }
        return lines;
    }

    // -------------------------------------------------------------------------
    // Temporary multi-child node for lenient parsing
    // -------------------------------------------------------------------------

    /**
     * Internal node used only during parsing — supports any number of children.
     * Converted to binary TreeNode after validation.
     */
    private static class RawNode {
        int taxonId = -1;                          // leaf: taxon ID; internal: -1
        final List<RawNode> children = new ArrayList<>();
        boolean isLeaf() { return children.isEmpty(); }
    }

    // -------------------------------------------------------------------------
    // Pass 2 – full parse
    // -------------------------------------------------------------------------

    /** Sentinel object pushed onto the stack to mark an open parenthesis. */
    private static final Object SENTINEL = new Object();

    /** Established unfiltered parser; kept free of allow-list branches. */
    private static Tree parseNewick(String s, int treeIdx, TaxonRegistry reg,
                                    int[] rootingCounts, boolean keepPolytomy) {
        int n = s.length();
        Deque<Object> stack = new ArrayDeque<>();   // contains RawNode or SENTINEL
        int i = 0;

        while (i < n) {
            char c = s.charAt(i);

            if (c == '(') {
                stack.push(SENTINEL);
                i++;

            } else if (c == ')') {
                // Collect all children pushed since the matching '('
                List<RawNode> children = new ArrayList<>();
                while (stack.peek() != SENTINEL) children.add((RawNode) stack.pop());
                stack.pop();   // remove sentinel

                // children were pushed left-to-right, popped right-to-left; restore order
                Collections.reverse(children);
                RawNode node = new RawNode();
                node.children.addAll(children);
                stack.push(node);

                i++;
                // skip optional internal label (e.g. bootstrap) then branch length
                i = skipLabelAndBranchLen(s, i, n);

            } else if (c == ',') {
                i++;

            } else if (c == ';') {
                break;

            } else if (c == ':') {
                // shouldn't appear at top level, but be safe
                i = skipBranchLen(s, i + 1, n);

            } else {
                // Leaf taxon name
                int start = i;
                i = scanLabelEnd(s, i, n);
                String name = decodeNewickLabel(s.substring(start, i));
                if (!name.isEmpty()) {
                    RawNode leaf = new RawNode();
                    leaf.taxonId = reg.getId(name);
                    stack.push(leaf);
                }
                i = skipBranchLen(s, i, n);   // skip ':length' if present
            }
        }

        if (stack.size() != 1) {
            throw new RuntimeException("Tree " + treeIdx
                + ": malformed Newick, stack size=" + stack.size());
        }
        RawNode rawRoot = (RawNode) stack.pop();
        if (rawRoot.isLeaf()) {
            throw new RuntimeException("Tree " + treeIdx + ": root is a leaf");
        }

        return buildTree(rawRoot, treeIdx, reg, rootingCounts, keepPolytomy);
    }

    /**
     * Allow-list parser. Unknown leaves disappear at tokenization time and
     * empty/unary clades are removed before the retained RawNode is converted.
     */
    private static Tree parseNewickRestricted(String s, int treeIdx,
                                              TaxonRegistry reg,
                                              int[] rootingCounts,
                                              boolean keepPolytomy) {
        int n = s.length();
        Deque<Object> stack = new ArrayDeque<>();
        int i = 0;

        while (i < n) {
            char c = s.charAt(i);
            if (c == '(') {
                stack.push(SENTINEL);
                i++;
            } else if (c == ')') {
                List<RawNode> children = new ArrayList<>();
                while (stack.peek() != SENTINEL) children.add((RawNode) stack.pop());
                stack.pop();
                Collections.reverse(children);

                if (children.size() == 1) {
                    stack.push(children.get(0));
                } else if (children.size() > 1) {
                    RawNode node = new RawNode();
                    node.children.addAll(children);
                    stack.push(node);
                }

                i++;
                i = skipLabelAndBranchLen(s, i, n);
            } else if (c == ',') {
                i++;
            } else if (c == ';') {
                break;
            } else if (c == ':') {
                i = skipBranchLen(s, i + 1, n);
            } else {
                int start = i;
                i = scanLabelEnd(s, i, n);
                String name = decodeNewickLabel(s.substring(start, i));
                int taxonId = name.isEmpty() ? -1 : reg.findId(name);
                if (taxonId >= 0) {
                    RawNode leaf = new RawNode();
                    leaf.taxonId = taxonId;
                    stack.push(leaf);
                }
                i = skipBranchLen(s, i, n);
            }
        }

        if (stack.isEmpty()) return null;
        if (stack.size() != 1) {
            throw new RuntimeException("Tree " + treeIdx
                + ": malformed Newick after taxon restriction, stack size=" + stack.size());
        }
        RawNode rawRoot = (RawNode) stack.pop();
        if (rawRoot.isLeaf()) return null;

        Tree tree = buildTree(rawRoot, treeIdx, reg, rootingCounts, keepPolytomy);
        boolean[] seen = new boolean[reg.size()];
        for (int taxonId : tree.postorderArray) {
            if (seen[taxonId]) {
                throw new IllegalArgumentException("Tree " + treeIdx
                    + " contains duplicate selected taxon: " + reg.getName(taxonId));
            }
            seen[taxonId] = true;
        }
        return tree;
    }

    private static Tree buildTree(RawNode rawRoot, int treeIdx, TaxonRegistry reg,
                                  int[] rootingCounts, boolean keepPolytomy) {
        int totalTaxa = reg.size();

        // Validate rooted arity and convert RawNode → TreeNode.
        int polytomyCountBefore = rootingCounts[2];
        TreeNode root = validateAndConvert(
            rawRoot, treeIdx, true, rootingCounts, keepPolytomy);
        boolean hasPolytomy = rootingCounts[2] != polytomyCountBefore;

        // Assign ranges and build postorderArray in one left-to-right DFS
        int[] postorderArray = new int[reg.size()]; // upper bound; trimmed below
        int[] counter = {0};
        assignRangesAndFillArray(root, postorderArray, counter);
        int leafCount = counter[0];
        postorderArray = Arrays.copyOf(postorderArray, leafCount);

        // Build inverse map
        int[] positionMap = new int[totalTaxa];
        Arrays.fill(positionMap, -1);
        for (int j = 0; j < leafCount; j++) positionMap[postorderArray[j]] = j;

        return new Tree(treeIdx, root, postorderArray, positionMap, leafCount, totalTaxa,
            hasPolytomy);
    }

    /** Exact child-order semantics of one clean.py deterministic first-pair refinement. */
    private static void resolveNodeFirstPair(RawNode node) {
        ArrayDeque<RawNode> work = new ArrayDeque<>(node.children);
        while (work.size() > 2) {
            RawNode first = work.removeFirst();
            RawNode second = work.removeFirst();
            RawNode joined = new RawNode();
            joined.children.add(first);
            joined.children.add(second);
            work.addLast(joined);
        }
        node.children.clear();
        node.children.addAll(work);
    }

    /**
     * Recursively validates a RawNode tree and converts it to a TreeNode:
     *
     *   isRoot=true,  2 children  → explicitly rooted binary tree.
     *   isRoot=true,  other arity → reject instead of inventing a root.
     *   isRoot=false, 2 children  → normal binary internal node.
     *   isRoot=false, ≥3 children → polytomous internal node (children[] array).
     *   leaf                      → leaf node.
     */
    private static TreeNode validateAndConvert(RawNode raw, int treeIdx, boolean isRoot,
                                               int[] rootingCounts, boolean keepPolytomy) {
        if (raw.isLeaf()) {
            TreeNode leaf = new TreeNode();
            leaf.taxonId = raw.taxonId;
            return leaf;
        }

        int nc = raw.children.size();
        if (isRoot && nc != 2) {
            throw new IllegalArgumentException("Tree " + treeIdx
                + " is not an explicitly rooted binary Newick tree: the top-level "
                + "node has " + nc + " children (STELAR-X requires exactly 2 and "
                + "never roots input trees arbitrarily)");
        }
        if (nc > (isRoot ? 3 : 2)) rootingCounts[4]++;
        if (isRoot && nc == 3) rootingCounts[5]++;
        if (!keepPolytomy && nc > 2) {
            rootingCounts[3]++;
            resolveNodeFirstPair(raw);
            nc = 2;
        }

        if (nc == 2) {
            TreeNode node = new TreeNode();
            node.left  = validateAndConvert(
                raw.children.get(0), treeIdx, false, rootingCounts, keepPolytomy);
            node.right = validateAndConvert(
                raw.children.get(1), treeIdx, false, rootingCounts, keepPolytomy);
            node.left.parent  = node;
            node.right.parent = node;
            return node;

        } else {
            // nc >= 3 && !isRoot  → polytomous internal node.
            rootingCounts[2]++;
            TreeNode node = new TreeNode();
            node.children = new TreeNode[nc];
            for (int j = 0; j < nc; j++) {
                TreeNode cj = validateAndConvert(
                    raw.children.get(j), treeIdx, false, rootingCounts, keepPolytomy);
                node.children[j] = cj;
                cj.parent = node;
            }
            node.left  = node.children[0];
            node.right = node.children[nc - 1];
            return node;
        }
    }

    /**
     * Single left-to-right DFS:
     *   - Leaf: assign rangeStart=counter, rangeEnd=counter+1, increment counter,
     *           write taxonId into postorderArray[counter].
     *   - Internal (binary or polytomous): recurse into children in order; range
     *     spans from leftmost (left) to rightmost (right) child.
     */
    private static void assignRangesAndFillArray(TreeNode node,
                                                  int[] arr, int[] counter) {
        if (node.isLeaf()) {
            node.rangeStart = counter[0];
            node.rangeEnd   = counter[0] + 1;
            arr[counter[0]] = node.taxonId;
            counter[0]++;
            return;
        }
        if (node.isPolytomous()) {
            for (TreeNode child : node.children) assignRangesAndFillArray(child, arr, counter);
        } else {
            assignRangesAndFillArray(node.left,  arr, counter);
            assignRangesAndFillArray(node.right, arr, counter);
        }
        node.rangeStart = node.left.rangeStart;
        node.rangeEnd   = node.right.rangeEnd;
    }

    private static void validateCompleteTaxonSet(Tree tree, TaxonRegistry reg, String inputFile) {
        int n = reg.size();
        if (tree.leafCount != n) {
            throw new IllegalArgumentException("Species tree taxon count (" + tree.leafCount
                + ") does not match gene-tree taxon count (" + n + "): " + inputFile);
        }

        boolean[] seen = new boolean[n];
        for (int taxon : tree.postorderArray) {
            if (seen[taxon]) {
                throw new IllegalArgumentException("Species tree contains duplicate taxon: "
                    + reg.getName(taxon));
            }
            seen[taxon] = true;
        }
        for (int i = 0; i < n; i++) {
            if (!seen[i]) {
                throw new IllegalArgumentException("Species tree is missing taxon: " + reg.getName(i));
            }
        }
    }

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /** True for characters that delimit a name or branch-length token. */
    private static boolean isDelim(char c) {
        return c == '(' || c == ')' || c == ',' || c == ':' || c == ';';
    }

    /**
     * Return the first structural delimiter after a possibly quoted Newick
     * label. Single quotes are the Newick standard; double quotes are accepted
     * as a compatibility extension. A doubled quote inside a quoted label is an
     * escaped literal quote.
     */
    private static int scanLabelEnd(String s, int i, int n) {
        char quote = 0;
        while (i < n) {
            char c = s.charAt(i);
            if (quote != 0) {
                if (c == quote) {
                    if (i + 1 < n && s.charAt(i + 1) == quote) {
                        i += 2;
                        continue;
                    }
                    quote = 0;
                }
                i++;
                continue;
            }
            if (c == '\'' || c == '"') {
                quote = c;
                i++;
                continue;
            }
            if (isDelim(c)) break;
            i++;
        }
        if (quote != 0) {
            throw new IllegalArgumentException("Unterminated quoted Newick label");
        }
        return i;
    }

    /**
     * Convert a Newick token to its logical taxon name. Thus {@code bird_name}
     * and {@code 'bird_name'} resolve to the same registry key, while doubled
     * quotes inside quoted names are unescaped.
     */
    private static String decodeNewickLabel(String raw) {
        String token = raw.trim();
        if (token.isEmpty()) return token;
        char quote = token.charAt(0);
        if (quote != '\'' && quote != '"') return token;

        StringBuilder decoded = new StringBuilder(token.length() - 2);
        for (int i = 1; i < token.length(); i++) {
            char c = token.charAt(i);
            if (c != quote) {
                decoded.append(c);
                continue;
            }
            if (i + 1 < token.length() && token.charAt(i + 1) == quote) {
                decoded.append(quote);
                i++;
                continue;
            }
            if (!token.substring(i + 1).trim().isEmpty()) {
                throw new IllegalArgumentException(
                    "Unexpected characters after quoted Newick label: " + token);
            }
            return decoded.toString();
        }
        throw new IllegalArgumentException("Unterminated quoted Newick label: " + token);
    }

    /** Skip digits/dots/e/+/- that make up a branch length value. */
    private static int skipBranchLen(String s, int i, int n) {
        if (i < n && s.charAt(i) == ':') i++;
        while (i < n && !isDelim(s.charAt(i))) i++;
        return i;
    }

    /**
     * After a ')' we may have: optional label (bootstrap or name), optional ':len'.
     * Skip both.
     */
    private static int skipLabelAndBranchLen(String s, int i, int n) {
        // Skip an optional label, respecting delimiters inside quoted names.
        i = scanLabelEnd(s, i, n);
        return skipBranchLen(s, i, n);
    }
}
