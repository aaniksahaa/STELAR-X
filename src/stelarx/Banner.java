package stelarx;

import stelarx.gpu.GPUWeightCalculator;

import java.io.PrintStream;

/**
 * Startup banner: system info + run configuration.
 * Written to stderr so it appears even when stdout is redirected.
 * ANSI colours are enabled only when the target stream is a real terminal
 * (or FORCE_COLOR is set), and suppressed when NO_COLOR is set.
 */
public class Banner {

    // ── ANSI colour codes ─────────────────────────────────────────────────────
    private static final String RST  = "\033[0m";
    private static final String BOLD = "\033[1m";
    private static final String DIM  = "\033[2m";
    private static final String CYAN = "\033[36m";
    private static final String GRN  = "\033[32m";
    private static final String YLW  = "\033[33m";
    private static final String WHT  = "\033[97m";
    private static final int TITLE_WIDTH = 63;

    private static final boolean USE_COLOR = detectColor(2);

    /** Exposed so PhaseLogger and other classes can share the same colour decision. */
    public static boolean useColor() { return USE_COLOR; }

    private static boolean detectColor(int fileDescriptor) {
        if (System.getenv("NO_COLOR")    != null) return false;
        if (System.getenv("FORCE_COLOR") != null) return true;
        // Check the actual target descriptor first so redirected output remains
        // machine-readable even when another stream is still attached to a terminal.
        try {
            String fd2 = java.nio.file.Files.readSymbolicLink(
                java.nio.file.Paths.get("/proc/self/fd/" + fileDescriptor)).toString();
            return fd2.startsWith("/dev/pts") || fd2.startsWith("/dev/tty");
        } catch (Exception e) {
            // /proc is unavailable on some platforms; use Java's console signal.
            return System.console() != null;
        }
    }

    private static String c(String code, String text) {
        return USE_COLOR ? code + text + RST : text;
    }

    private static String c(boolean color, String code, String text) {
        return color ? code + text + RST : text;
    }

    private static String fmtMiB(long mib) {
        return mib >= 1024 ? String.format("%.1f GB", mib / 1024.0) : mib + " MB";
    }

    // ── Public entry point ────────────────────────────────────────────────────

    /** Print the shared title box, using the colour state of the target stream. */
    public static void printTitle(PrintStream out) {
        boolean color = out == System.out ? detectColor(1) : USE_COLOR;
        printTitle(out, color);
    }

    /** Print the title box and compact version greeting to stdout. */
    public static void printVersion() {
        boolean color = detectColor(1);
        printTitle(System.out, color);
        System.out.println(c(color, WHT,
            "Welcome to STELAR-X version " + Main.VERSION + "!"));
        System.out.println();
    }

    private static void printTitle(PrintStream out, boolean color) {
        // w = number of ═ characters (= total inner + 2 for the space padding each side)
        String title = "STELAR-X  v" + Main.VERSION;
        int inner = TITLE_WIDTH - 2;
        int lpad  = (inner - title.length()) / 2;
        int rpad  = inner - title.length() - lpad;
        String paddedTitle = " ".repeat(Math.max(0, lpad))
                           + title
                           + " ".repeat(Math.max(0, rpad));

        out.println();
        out.println("  " + c(color, BOLD + CYAN, "╔" + "═".repeat(TITLE_WIDTH) + "╗"));
        out.println("  " + c(color, BOLD + CYAN, "║") + " "
                         + c(color, BOLD + WHT, paddedTitle)
                         + " " + c(color, BOLD + CYAN, "║"));
        out.println("  " + c(color, BOLD + CYAN, "╚" + "═".repeat(TITLE_WIDTH) + "╝"));
        out.println();
    }

    public static void print(Config cfg) {
        PrintStream out = System.err;
        printTitle(out, USE_COLOR);

        String sep = "─".repeat(TITLE_WIDTH - 2);

        // ── System section ─────────────────────────────────────────────────
        out.println("  " + c(BOLD, "System") + "  " + c(DIM, sep.substring(0, sep.length() - 4)));

        int available = Runtime.getRuntime().availableProcessors();
        int using     = cfg.getThreadCount();
        out.println("    " + String.format("%-8s %s  →  %s",
            "CPU",
            c(WHT, available + " cores available"),
            c(available == using ? WHT : YLW, using + " threads configured")));

        GPUWeightCalculator.Probe gpuProbe = GPUWeightCalculator.probe();
        if (gpuProbe.cudaAvailable()) {
            String device = gpuProbe.deviceName() + "  (CC "
                + gpuProbe.computeMajor() + "." + gpuProbe.computeMinor() + ")";
            out.println("    " + String.format("%-8s %s  ·  %s total  ·  %s free%s",
                "GPU",
                c(WHT, device),
                c(WHT, fmtMiB(gpuProbe.totalMiB())),
                c(gpuProbe.freeMiB() > gpuProbe.totalMiB() / 4 ? GRN : YLW,
                    fmtMiB(gpuProbe.freeMiB())),
                "  " + c(GRN, "✓ CUDA usable")));
        } else {
            out.println("    " + String.format("%-8s %s",
                "GPU", c(YLW, "unavailable  (CPU fallback ready)")));
            if (cfg.getRequestedComputeMode() != Config.ComputeMode.CPU) {
                out.println("    " + String.format("%-8s %s", "",
                    c(DIM, gpuProbe.detail())));
            }
        }
        out.println();

        // ── Run configuration section ───────────────────────────────────────
        out.println("  " + c(BOLD, "Run Configuration") + "  " + c(DIM, sep.substring(0, sep.length() - 15)));
        out.println();

        boolean gpuMode = cfg.getComputeMode() == Config.ComputeMode.GPU;
        String naTag    = c(DIM, "(n/a)");

        // ── I/O ────────────────────────────────────────────────────────────
        String inputPath = cfg.getInputFile();
        String displayInput = "(none)";
        if (inputPath != null) {
            String[] parts = inputPath.replace('\\', '/').split("/");
            displayInput = parts.length > 2
                ? "…/" + parts[parts.length - 2] + "/" + parts[parts.length - 1]
                : inputPath;
        }
        out.println("    " + row("Input file",    c(WHT, displayInput)));
        String outputPath = cfg.getOutputFile();
        out.println("    " + row("Output file",   outputPath != null
                                                  ? c(WHT, outputPath)
                                                  : c(DIM, "(stdout)")));
        String logPath = cfg.getLogFile();
        if (logPath != null) {
            out.println("    " + row("Terminal log", c(WHT, logPath)));
        }
        if (cfg.getTaxaFile() != null) {
            out.println("    " + row("Taxon allow-list", c(WHT, cfg.getTaxaFile())));
        }
        out.println();

        // ── Compute ────────────────────────────────────────────────────────
        String computeValue = gpuMode ? c(GRN, "GPU") : c(WHT, "CPU");
        computeValue += c(DIM, "  (" + cfg.getComputeModeDetail() + ")");
        out.println("    " + row("Compute mode", computeValue));
        out.println("    " + row("CPU threads",    c(available == using ? WHT : YLW, String.valueOf(using))
                                                  + c(DIM, "  (" + available + " available)")));
        out.println("    " + row("Tree treatment", c(WHT, cfg.getTreatAsUnrooted() ? "unrooted" : "rooted")));
        String inferencePolytomies = cfg.isScoreOnly()
            ? c(DIM, "(n/a; no inference)")
            : c(WHT, cfg.isKeepPolytomyDuringInference()
                ? "keep"
                : "resolve  (deterministic first-pair refinement)");
        out.println("    " + row("Inference polytomies", inferencePolytomies));
        out.println("    " + row("Triplet-score polytomies",
            c(WHT, "keep  (always; native input topology)")));
        out.println();

        // ── Search ─────────────────────────────────────────────────────────
        out.println("    " + row("Search mode",   c(WHT, cfg.getSearchMode().name().toLowerCase())));
        out.println("    " + row("Hash seeds",    c(WHT, String.valueOf(cfg.getNumHashSeeds()))));

        String verbStr = switch (cfg.getVerbosity()) {
            case Logging.QUIET -> c(DIM, "quiet");
            case Logging.DEBUG -> c(YLW, "debug");
            case Logging.TRACE -> c(YLW, "trace");
            default            -> c(WHT, "info");
        };
        out.println("    " + row("Verbosity",     verbStr));
        out.println();

        // ── GPU parameters ─────────────────────────────────────────────────
        // Weight batching mode
        String batchStr;
        if (!gpuMode) {
            batchStr = naTag;
        } else if (!cfg.isGpuBatch()) {
            batchStr = c(WHT, "disabled  (single launch)");
        } else if (cfg.getGpuNumBatches() > 0) {
            batchStr = c(WHT, cfg.getGpuNumBatches() + " batches") + c(DIM, "  (manual)");
        } else if (cfg.getGpuBatchSize() > 0) {
            batchStr = c(WHT, "batch-size " + cfg.getGpuBatchSize()) + c(DIM, "  (manual)");
        } else if (cfg.isGpuVramControlFactorSet()) {
            batchStr = c(WHT, "vram-control-factor") + c(DIM, "  (resident-relative)");
        } else {
            batchStr = c(WHT, "auto") + c(DIM, "  (method-specific)");
        }
        out.println("    " + row("Weight batching",          batchStr));

        // Batching sub-parameter: show the active knob
        if (gpuMode && cfg.isGpuBatch()) {
            if (cfg.isGpuVramControlFactorSet()) {
                out.println("    " + row("  VRAM control factor",
                    c(WHT, String.format("%.3f", cfg.getGpuVramControlFactor()))
                    + c(DIM, "  →  mem(batch) = F × mem(resident)")));
            } else if (cfg.getGpuNumBatches() <= 0 && cfg.getGpuBatchSize() <= 0) {
                out.println("    " + row("  VRAM occupancy factor",
                    c(WHT, String.format("%.2f", cfg.getGpuVramFraction()))
                    + c(DIM, "  →  batch = freeVRAM × factor")));
                out.println("    " + row("  Tree-walk scratch cap",
                    c(WHT, cfg.getGpuTreeWalkVramCapMiB() + " MiB")
                    + c(DIM, "  (automatic simple-tree-walk batches)")));
            }
        }


        // DP state-space cap (only meaningful for GPU + FULL search)
        boolean dpRelevant = gpuMode && cfg.getSearchMode() == Config.SearchMode.FULL;
        out.println("    " + row("DP state-space construction memory cap (GPU)", dpRelevant
                ? c(WHT, fmtBytes(cfg.getGpuDpOutputCapBytes()))
                  + c(DIM, "  (" + cfg.getGpuDpOutputCapTriples() + " triples)")
                : naTag));

        out.println();
        out.println("  " + c(DIM, "─".repeat(TITLE_WIDTH)));
        out.println();
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    // Label column width for Run Configuration rows
    private static String row(String label, String value) {
        return String.format("%-46s %s", c(DIM, label), value);
    }

    private static String fmtBytes(long bytes) {
        if (bytes >= 1_000_000_000L) return String.format("%.1f GB", bytes / 1e9);
        if (bytes >= 1_000_000L)     return String.format("%.0f MB", bytes / 1e6);
        if (bytes >= 1_000L)         return String.format("%.0f KB", bytes / 1e3);
        return bytes + " B";
    }
}
