package utils;

public class Logging {
    public static final int QUIET = 0;
    public static final int INFO  = 1;
    public static final int DEBUG = 2;
    public static final int TRACE = 3;

    private static int level = INFO;
    private static final long startTime = System.nanoTime();

    public static void setLevel(int l) { level = l; }
    public static int  getLevel()      { return level; }
    public static boolean isDebug()    { return level >= DEBUG; }
    public static boolean isTrace()    { return level >= TRACE; }

    public static void info (String fmt, Object... a) { if (level>=INFO)  log("INFO ", fmt, a); }
    public static void debug(String fmt, Object... a) { if (level>=DEBUG) log("DEBUG", fmt, a); }
    public static void trace(String fmt, Object... a) { if (level>=TRACE) log("TRACE", fmt, a); }
    public static void error(String fmt, Object... a) { log("ERROR", fmt, a); }

    private static void log(String tag, String fmt, Object[] a) {
        long ms = (System.nanoTime() - startTime) / 1_000_000;
        System.err.printf("[%5dms][%s] %s%n", ms, tag, String.format(fmt, a));
    }
}
