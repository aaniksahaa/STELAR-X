package stelarx;

import java.io.IOException;
import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

/**
 * Duplicates Java stdout and stderr into one run log, excluding transient
 * carriage-return progress repaints.
 *
 * The portable Unix launcher performs the duplication outside the JVM so that
 * output written directly by CUDA/native libraries is captured as well.  This
 * class is the fallback for direct Java and Windows-launcher invocations.
 */
final class TerminalLog implements AutoCloseable {
    static final class SetupException extends Exception {
        SetupException(String message) { super(message); }
        SetupException(String message, Throwable cause) { super(message, cause); }
    }

    private final PrintStream originalOut;
    private final PrintStream originalErr;
    private final OutputStream log;
    private boolean closed;

    private TerminalLog(PrintStream originalOut, PrintStream originalErr,
                        OutputStream log) {
        this.originalOut = originalOut;
        this.originalErr = originalErr;
        this.log = log;
    }

    static TerminalLog installFromArgs(String[] args) throws SetupException {
        if ("1".equals(System.getenv("STELARX_LOG_CAPTURED"))) return null;

        String logName = optionValue(args, "--log-file");
        if (logName == null) return null;
        if (logName.isBlank()) {
            throw new SetupException("--log-file requires a non-empty file path.");
        }

        Path logPath = Path.of(logName).toAbsolutePath().normalize();
        protectAnalysisFiles(args, logPath);

        try {
            Path parent = logPath.getParent();
            if (parent != null) Files.createDirectories(parent);
            OutputStream file = Files.newOutputStream(logPath,
                StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING,
                StandardOpenOption.WRITE);
            OutputStream sink = new ProgressFilteringOutputStream(file);
            PrintStream oldOut = System.out;
            PrintStream oldErr = System.err;
            TerminalLog session = new TerminalLog(oldOut, oldErr, sink);
            System.setOut(new PrintStream(new TeeOutputStream(oldOut, sink), true,
                StandardCharsets.UTF_8));
            System.setErr(new PrintStream(new TeeOutputStream(oldErr, sink), true,
                StandardCharsets.UTF_8));
            Runtime.getRuntime().addShutdownHook(
                new Thread(session::closeQuietly, "stelarx-log-closer"));
            return session;
        } catch (IOException | SecurityException e) {
            throw new SetupException("cannot open log file '" + logPath + "': "
                + e.getMessage(), e);
        }
    }

    private static String optionValue(String[] args, String option)
            throws SetupException {
        String value = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(option)) {
                if (++i >= args.length) {
                    throw new SetupException(option + " requires a file path.");
                }
                value = args[i];
            } else if (args[i].startsWith(option + "=")) {
                value = args[i].substring(option.length() + 1);
            }
        }
        return value;
    }

    private static void protectAnalysisFiles(String[] args, Path logPath)
            throws SetupException {
        for (int i = 0; i < args.length; i++) {
            String option = args[i];
            if ((option.equals("-i") || option.equals("--input")
                    || option.equals("-o") || option.equals("--output")
                    || option.equals("--score-species-tree")
                    || option.equals("--taxa-file")
                    || option.equals("--species-list")
                    || option.equals("--species-list-file")
                    || option.equals("--dump-clusters")
                    || option.equals("--dump-completed-gene-trees"))
                    && i + 1 < args.length) {
                Path dataPath = Path.of(args[++i]).toAbsolutePath().normalize();
                if (logPath.equals(dataPath)) {
                    throw new SetupException("--log-file must differ from every analysis input and output: "
                        + logPath);
                }
            }
        }
    }

    @Override
    public synchronized void close() {
        if (closed) return;
        closed = true;
        System.out.flush();
        System.err.flush();
        System.setOut(originalOut);
        System.setErr(originalErr);
        try {
            log.close();
        } catch (IOException ignored) {
            // The run has already finished; there is nowhere safer to report this.
        }
    }

    private void closeQuietly() {
        try { close(); } catch (Throwable ignored) { }
    }

    private static final class TeeOutputStream extends OutputStream {
        private final OutputStream terminal;
        private final OutputStream log;

        TeeOutputStream(OutputStream terminal, OutputStream log) {
            this.terminal = terminal;
            this.log = log;
        }

        @Override
        public void write(int value) throws IOException {
            synchronized (log) {
                terminal.write(value);
                log.write(value);
            }
        }

        @Override
        public void write(byte[] bytes, int offset, int length) throws IOException {
            synchronized (log) {
                terminal.write(bytes, offset, length);
                log.write(bytes, offset, length);
            }
        }

        @Override
        public void flush() throws IOException {
            synchronized (log) {
                terminal.flush();
                log.flush();
            }
        }
    }

    /** Omits carriage-return repaint records while retaining normal log lines. */
    private static final class ProgressFilteringOutputStream extends OutputStream {
        private static final byte[] PROGRESS_GLYPH = "▸".getBytes(StandardCharsets.UTF_8);

        private final OutputStream file;
        private final ByteArrayOutputStream line = new ByteArrayOutputStream(256);
        private boolean hasCarriageReturn;

        ProgressFilteringOutputStream(OutputStream file) {
            this.file = file;
        }

        @Override
        public void write(int value) throws IOException {
            int b = value & 0xff;
            if (b == '\n') {
                finishLine(true);
                return;
            }
            if (b == '\r') hasCarriageReturn = true;
            line.write(b);
        }

        @Override
        public void write(byte[] bytes, int offset, int length) throws IOException {
            for (int i = offset; i < offset + length; i++) write(bytes[i]);
        }

        private void finishLine(boolean newline) throws IOException {
            byte[] bytes = line.toByteArray();
            boolean progress = hasCarriageReturn
                && (!isSingleTrailingCarriageReturn(bytes) || contains(bytes, PROGRESS_GLYPH));
            if (!progress) {
                int length = bytes.length;
                if (length > 0 && bytes[length - 1] == '\r') length--;
                file.write(bytes, 0, length);
                if (newline) file.write('\n');
            }
            line.reset();
            hasCarriageReturn = false;
        }

        private static boolean isSingleTrailingCarriageReturn(byte[] bytes) {
            int count = 0;
            for (byte b : bytes) if (b == '\r') count++;
            return count == 1 && bytes.length > 0 && bytes[bytes.length - 1] == '\r';
        }

        private static boolean contains(byte[] bytes, byte[] needle) {
            outer:
            for (int i = 0; i <= bytes.length - needle.length; i++) {
                for (int j = 0; j < needle.length; j++) {
                    if (bytes[i + j] != needle[j]) continue outer;
                }
                return true;
            }
            return false;
        }

        @Override
        public void flush() throws IOException {
            file.flush();
        }

        @Override
        public void close() throws IOException {
            if (line.size() > 0) finishLine(false);
            file.close();
        }
    }
}
