package stelarx;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Stream;

/** Verifies that fatal reports are isolated in a dedicated directory. */
public final class FatalReporterTest {
    public static void main(String[] args) throws Exception {
        if (args.length != 1) throw new IllegalArgumentException("pass a temporary directory");
        Path crashDir = Path.of(args[0]).resolve("nested").resolve("crash_logs");
        System.setProperty(FatalReporter.CRASH_DIR_PROPERTY, crashDir.toString());

        FatalReporter.report(new IllegalStateException("crash-directory-sentinel"),
            new String[] {"--diagnose"});

        check(Files.isDirectory(crashDir), "crash directory was not created");
        List<Path> reports;
        try (Stream<Path> files = Files.list(crashDir)) {
            reports = files.filter(path -> path.getFileName().toString()
                .matches("stelarx-crash-[0-9]{8}-[0-9]{6}-[0-9]+[.]log"))
                .toList();
        }
        check(reports.size() == 1, "expected one crash report, found " + reports);
        String text = Files.readString(reports.get(0));
        check(text.contains("STELAR-X FATAL ERROR"), "missing report header");
        check(text.contains("crash-directory-sentinel"), "missing failure message");
        check(text.contains("Command: stelarx --diagnose"), "missing command line");
        System.out.println("Crash-report directory isolation: PASS");
    }

    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }
}
