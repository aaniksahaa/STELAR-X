package stelarx;

import stelarx.gpu.GPUWeightCalculator;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;

/** Best-effort fatal report for Java exceptions and memory failures. */
public final class FatalReporter {
    static final String CRASH_DIR_PROPERTY = "stelarx.crashDir";
    static final String CRASH_DIR_ENV = "STELARX_CRASH_DIR";

    private FatalReporter() {}

    /** Best-effort early creation, primarily for packaged JVM fatal-error logs. */
    public static void prepareCrashDirectory() {
        try {
            createCrashDirectory();
        } catch (Throwable ignored) {
            // A real failure will retry and print the concrete write error.
        }
    }

    public static void report(Throwable failure, String[] args) {
        String text = build(failure, args);
        System.err.println(text);

        String stamp = DateTimeFormatter.ofPattern("yyyyMMdd-HHmmss")
            .withZone(ZoneOffset.UTC).format(Instant.now());
        String name = "stelarx-crash-" + stamp + "-" + ProcessHandle.current().pid() + ".log";
        try {
            Path report = createCrashDirectory().resolve(name);
            Files.writeString(report, text + System.lineSeparator(), StandardCharsets.UTF_8);
            System.err.println("Crash report written to: " + report.toAbsolutePath());
        } catch (Throwable writeFailure) {
            System.err.println("Could not write crash report: " + writeFailure.getMessage());
        }
    }

    /**
     * Creates the configured crash directory on demand. The default is a
     * {@code crash_logs} child of the process working directory. Launchers can
     * override it with {@code -Dstelarx.crashDir} or {@code STELARX_CRASH_DIR}.
     */
    static Path createCrashDirectory() throws IOException {
        String configured = System.getProperty(CRASH_DIR_PROPERTY);
        if (configured == null || configured.isBlank()) {
            configured = System.getenv(CRASH_DIR_ENV);
        }
        Path directory = configured == null || configured.isBlank()
            ? Path.of(System.getProperty("user.dir", "."), "crash_logs")
            : Path.of(configured);
        directory = directory.toAbsolutePath().normalize();
        try {
            return ensureDirectory(directory);
        } catch (IOException primaryFailure) {
            Path fallback = Path.of(System.getProperty("java.io.tmpdir", "."),
                "stelarx-crash_logs").toAbsolutePath().normalize();
            if (fallback.equals(directory)) throw primaryFailure;
            try {
                return ensureDirectory(fallback);
            } catch (IOException fallbackFailure) {
                fallbackFailure.addSuppressed(primaryFailure);
                throw fallbackFailure;
            }
        }
    }

    private static Path ensureDirectory(Path directory) throws IOException {
        Files.createDirectories(directory);
        if (!Files.isDirectory(directory)) {
            throw new IOException("crash-log path is not a directory: " + directory);
        }
        return directory;
    }

    private static String build(Throwable failure, String[] args) {
        Runtime rt = Runtime.getRuntime();
        StringWriter sw = new StringWriter();
        PrintWriter out = new PrintWriter(sw);
        out.println();
        out.println("================ STELAR-X FATAL ERROR ================");
        out.println("Version: " + Main.VERSION);
        out.println("Time (UTC): " + Instant.now());
        out.println("Phase: " + PhaseLogger.currentPhase());
        out.println("Failure: " + failure.getClass().getName());
        out.println("Message: " + String.valueOf(failure.getMessage()));
        out.println("OS/arch: " + System.getProperty("os.name") + " "
            + System.getProperty("os.version") + " / " + System.getProperty("os.arch"));
        out.println("Runtime: " + System.getProperty("java.runtime.version"));
        out.println("Memory: used=" + mib(rt.totalMemory() - rt.freeMemory())
            + " MiB, committed=" + mib(rt.totalMemory()) + " MiB, max=" + mib(rt.maxMemory()) + " MiB");
        GPUWeightCalculator.Probe gpu = GPUWeightCalculator.probe();
        out.println("CUDA: " + (gpu.cudaAvailable()
            ? gpu.deviceName() + " (CC " + gpu.computeMajor() + "." + gpu.computeMinor() + ")"
            : "unavailable: " + gpu.detail()));
        out.println("Command: stelarx " + quoteArgs(args));
        if (failure instanceof OutOfMemoryError) {
            out.println();
            out.println("Likely remedy: close other memory-heavy jobs, use a machine with more RAM,");
            out.println("or reduce search-space enrichment. The portable launcher allows the JVM");
            out.println("to use up to 85% of physical/container memory by default.");
        }
        out.println();
        failure.printStackTrace(out);
        out.println("======================================================");
        out.flush();
        return sw.toString();
    }

    private static long mib(long bytes) { return bytes / (1024L * 1024L); }

    private static String quoteArgs(String[] args) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            if (i > 0) sb.append(' ');
            String a = args[i];
            if (a.matches("[A-Za-z0-9_./:=+,-]+")) sb.append(a);
            else sb.append('\'').append(a.replace("'", "'\\''")).append('\'');
        }
        return sb.toString();
    }
}
