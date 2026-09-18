package stelarx.tree;

import stelarx.Config;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/** Parser-backed taxon-set extraction, allow-list loading, and coverage statistics. */
public final class TreeTaxa {
    private TreeTaxa() {}

    /**
     * One streaming scan of a one-Newick-per-line file.  Only the union and
     * intersection sets are retained; per-tree sets are released immediately.
     */
    public static Scan scan(String inputFile, Set<String> reference) throws IOException {
        LinkedHashSet<String> union = new LinkedHashSet<>();
        LinkedHashSet<String> intersection = null;
        long totalMissing = 0L;
        int minMissing = Integer.MAX_VALUE;
        int maxMissing = 0;
        int treeCount = 0;

        try (BufferedReader reader = Files.newBufferedReader(
                Path.of(inputFile), StandardCharsets.UTF_8)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                LinkedHashSet<String> current = new LinkedHashSet<>();
                TreeParser.forEachTaxonName(line, current::add);
                if (current.isEmpty()) {
                    throw new IllegalArgumentException(
                        "Tree " + treeCount + " contains no taxon names: " + inputFile);
                }

                union.addAll(current);
                if (intersection == null) intersection = new LinkedHashSet<>(current);
                else intersection.retainAll(current);

                if (reference != null) {
                    int present = 0;
                    for (String name : current) if (reference.contains(name)) present++;
                    int missing = reference.size() - present;
                    totalMissing += missing;
                    minMissing = Math.min(minMissing, missing);
                    maxMissing = Math.max(maxMissing, missing);
                }
                treeCount++;
            }
        }

        if (treeCount == 0) {
            throw new IllegalArgumentException("Tree file is empty: " + inputFile);
        }
        if (reference == null) {
            minMissing = 0;
            maxMissing = 0;
        }
        return new Scan(treeCount, union,
            intersection == null ? new LinkedHashSet<>() : intersection,
            totalMissing, minMissing, maxMissing);
    }

    /**
     * Coverage scan specialized for inference allow-lists. Only selected names
     * are retained, so memory is O(size of the allow-list) even when the input
     * contains a much larger outside taxon universe.
     */
    public static SelectionScan scanSelection(String inputFile, Set<String> selected)
            throws IOException {
        LinkedHashSet<String> selectedUnion = new LinkedHashSet<>();
        long totalMissing = 0L;
        long ignoredLeafOccurrences = 0L;
        int minMissing = Integer.MAX_VALUE;
        int maxMissing = 0;
        int treeCount = 0;

        try (BufferedReader reader = Files.newBufferedReader(
                Path.of(inputFile), StandardCharsets.UTF_8)) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;

                LinkedHashSet<String> current = new LinkedHashSet<>();
                long[] ignoredInTree = {0L};
                int[] leafCount = {0};
                TreeParser.forEachTaxonName(line, name -> {
                    leafCount[0]++;
                    if (selected.contains(name)) current.add(name);
                    else ignoredInTree[0]++;
                });
                if (leafCount[0] == 0) {
                    throw new IllegalArgumentException(
                        "Tree " + treeCount + " contains no taxon names: " + inputFile);
                }

                selectedUnion.addAll(current);
                int missing = selected.size() - current.size();
                totalMissing += missing;
                ignoredLeafOccurrences += ignoredInTree[0];
                minMissing = Math.min(minMissing, missing);
                maxMissing = Math.max(maxMissing, missing);
                treeCount++;
            }
        }

        if (treeCount == 0) {
            throw new IllegalArgumentException("Tree file is empty: " + inputFile);
        }
        return new SelectionScan(treeCount, selectedUnion, totalMissing,
            minMissing, maxMissing, ignoredLeafOccurrences);
    }

    /** Read a taxon allow-list: one name per non-empty line, retaining file order. */
    public static TaxaList readTaxaList(String taxaFile) throws IOException {
        LinkedHashSet<String> names = new LinkedHashSet<>();
        int nonEmptyLines = 0;
        try (BufferedReader reader = Files.newBufferedReader(
                Path.of(taxaFile), StandardCharsets.UTF_8)) {
            String line;
            while ((line = reader.readLine()) != null) {
                String name = line.trim();
                if (name.isEmpty()) continue;
                nonEmptyLines++;
                names.add(name);
            }
        }
        if (names.isEmpty()) {
            throw new IllegalArgumentException("Taxa file is empty: " + taxaFile);
        }
        return new TaxaList(names, nonEmptyLines - names.size());
    }

    /** Write a deterministic, delimiter-free taxa list (one sorted name per line). */
    public static int writeExtracted(String inputFile, String outputFile,
                                     Config.TaxaSetMode mode) throws IOException {
        Scan scan = scan(inputFile, null);
        List<String> names = new ArrayList<>(mode == Config.TaxaSetMode.INTERSECTION
            ? scan.intersection() : scan.union());
        Collections.sort(names);

        if (outputFile == null) {
            PrintWriter out = new PrintWriter(System.out, true, StandardCharsets.UTF_8);
            for (String name : names) out.println(name);
            out.flush();
        } else {
            Path output = Path.of(outputFile);
            Path parent = output.toAbsolutePath().getParent();
            if (parent != null) Files.createDirectories(parent);
            try (BufferedWriter writer = Files.newBufferedWriter(
                    output, StandardCharsets.UTF_8)) {
                for (String name : names) {
                    writer.write(name);
                    writer.newLine();
                }
            }
        }
        return names.size();
    }

    public record Scan(int treeCount,
                       LinkedHashSet<String> union,
                       LinkedHashSet<String> intersection,
                       long totalMissing,
                       int minMissing,
                       int maxMissing) {
        public double meanMissing() {
            return treeCount == 0 ? 0.0 : (double) totalMissing / treeCount;
        }
    }

    public record SelectionScan(int treeCount,
                                LinkedHashSet<String> selectedUnion,
                                long totalMissing,
                                int minMissing,
                                int maxMissing,
                                long ignoredLeafOccurrences) {
        public double meanMissing() {
            return treeCount == 0 ? 0.0 : (double) totalMissing / treeCount;
        }
    }

    public record TaxaList(LinkedHashSet<String> names, int duplicateLines) {}
}
