import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

import taxon.Taxon;
import tree.Tree;

public class CompareTrees {
    public static void main(String[] args) throws FileNotFoundException {
        if (args.length != 2) {
            System.out.println("Usage: java CompareTrees <tree_file_1> <tree_file_2>");
            System.exit(1);
        }
        String file1 = args[0];
        String file2 = args[1];

        String newick1 = readFirstNonEmptyLine(file1);
        String newick2 = readFirstNonEmptyLine(file2);

        if (newick1 == null || newick2 == null) {
            System.err.println("Error: One or both files are empty or missing a valid Newick tree.");
            System.exit(2);
        }

        // Build taxa map from both trees (union of taxa)
        Map<String, Taxon> taxaMap = new HashMap<>();
        for (String label : extractTaxa(newick1)) {
            taxaMap.put(label, new Taxon(label));
        }
        for (String label : extractTaxa(newick2)) {
            taxaMap.putIfAbsent(label, new Taxon(label));
        }

        Tree t1 = new Tree(newick1, taxaMap);
        Tree t2 = new Tree(newick2, taxaMap);

        boolean identical = t1.equalsStructure(t2);
        if (identical) {
            System.out.println("The trees are structurally identical (rooted, ignoring child order).");
        } else {
            System.out.println("The trees are NOT structurally identical.");
        }
    }

    private static String readFirstNonEmptyLine(String filePath) throws FileNotFoundException {
        Scanner scanner = new Scanner(new File(filePath));
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine().trim();
            if (!line.isEmpty()) {
                scanner.close();
                return line;
            }
        }
        scanner.close();
        return null;
    }

    private static java.util.Set<String> extractTaxa(String newick) {
        java.util.Set<String> taxa = new java.util.HashSet<>();
        int n = newick.length();
        int i = 0;
        while (i < n) {
            char c = newick.charAt(i);
            if (c == '(' || c == ')' || c == ',' || c == ';') {
                i++;
            } else {
                StringBuilder sb = new StringBuilder();
                while (i < n) {
                    char cc = newick.charAt(i);
                    if (cc == '(' || cc == ')' || cc == ',' || cc == ':' || cc == ';') break;
                    sb.append(cc);
                    i++;
                }
                String label = sb.toString();
                if (!label.isEmpty()) taxa.add(label);
                // skip branch length
                if (i < n && newick.charAt(i) == ':') {
                    i++;
                    while (i < n) {
                        char cc = newick.charAt(i);
                        if (cc == ',' || cc == ')' || cc == ';') break;
                        i++;
                    }
                }
            }
        }
        return taxa;
    }
} 