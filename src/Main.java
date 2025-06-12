import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.util.Scanner;
import utils.Config;

/**
 * Main entry point for basic phylogeny project setup.
 * 
 * This is a minimal implementation that reads gene trees from an input file
 * and outputs basic statistics to demonstrate the setup is working.
 */
public class Main {

    /**
     * Main method that processes gene trees and outputs basic statistics.
     * 
     * Command line arguments:
     * -i <input_file>  - Input file path containing gene trees in Newick format
     * -o <output_file> - Output file path for results
     */
    public static void main(String[] args) throws IOException {

        String inputFilePath = null;
        String outputFilePath = null;

        // Parse command line arguments
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i") && i + 1 < args.length) {
                inputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            } else if (args[i].equals("-o") && i + 1 < args.length) {
                outputFilePath = args[i + 1];
                i++; // Skip next argument as it's the file path
            }
        }

        // Validate required arguments
        if (inputFilePath == null || outputFilePath == null) {
            System.out.println("Usage: java Main -i <input_file> -o <output_file>");
            System.exit(-1);
        }

        System.out.println("Input file: " + inputFilePath);
        System.out.println("Output file: " + outputFilePath);

        try {
            // Count gene trees by reading lines from input file
            int geneTreeCount = 0;
            
            File inputFile = new File(inputFilePath);
            if (inputFile.exists()) {
                Scanner scanner = new Scanner(inputFile);
                while (scanner.hasNextLine()) {
                    String line = scanner.nextLine().trim();
                    if (line.length() > 0 && !line.startsWith("#")) {
                        geneTreeCount++;
                    }
                }
                scanner.close();
            } else {
                System.out.println("Input file does not exist, using default count of 0");
            }

            // Create output with spacing from Config
            StringBuilder output = new StringBuilder();
            
            // Add n spaces before the output
            for (int i = 0; i < Config.n; i++) {
                output.append(" ");
            }
            
            // Add the number of gene trees
            output.append("Number of gene trees: ").append(geneTreeCount);

            // Write to output file
            FileWriter writer = new FileWriter(outputFilePath);
            writer.write(output.toString());
            writer.close();

            System.out.println("Successfully processed " + geneTreeCount + " gene trees");
            System.out.println("Results written to: " + outputFilePath);

        } catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
            e.printStackTrace();
            System.exit(-1);
        }
    }
}
