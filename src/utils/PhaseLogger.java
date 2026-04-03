package utils;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Step-by-step phase logging with per-phase peak CPU RAM.
 *
 * A background thread runs for every phase:
 *   - cpu-ram-poller : reads /proc/self/status VmRSS every 50 ms (total process RSS,
 *                      includes JNI native memory — matches what htop shows)
 *
 * Output format:
 *   Phase N  Name  [CPU]
 *   ...output interleaved...
 *        1234 ms  |  peak RAM: 8192 MiB
 */
public class PhaseLogger {

    // ── CPU RAM poller (/proc/self/status VmRSS) ──────────────────────────────
    private static final AtomicLong peakCpuRssMiB   = new AtomicLong(0);
    private static volatile boolean cpuPolling       = false;
    private static Thread           cpuPollThread    = null;

    /** Read resident set size in MiB from /proc/self/status, or -1 on failure. */
    private static long readVmRssMiB() {
        try (java.io.BufferedReader br = new java.io.BufferedReader(
                new java.io.FileReader("/proc/self/status"))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("VmRSS:")) {
                    // format: "VmRSS:   12345 kB"
                    String[] parts = line.trim().split("\\s+");
                    if (parts.length >= 2) return Long.parseLong(parts[1]) / 1024;
                }
            }
        } catch (Exception ignored) {}
        // fallback: Java heap only
        Runtime rt = Runtime.getRuntime();
        return (rt.totalMemory() - rt.freeMemory()) / (1024 * 1024);
    }

    private static void startCpuPolling() {
        peakCpuRssMiB.set(0);
        cpuPolling = true;
        cpuPollThread = new Thread(() -> {
            while (cpuPolling) {
                long rss = readVmRssMiB();
                if (rss > 0) {
                    long cur;
                    do { cur = peakCpuRssMiB.get(); }
                    while (rss > cur && !peakCpuRssMiB.compareAndSet(cur, rss));
                }
                try { Thread.sleep(50); } catch (InterruptedException e) { break; }
            }
        }, "cpu-ram-poller");
        cpuPollThread.setDaemon(true);
        cpuPollThread.start();
    }

    /** Stop CPU polling and return a formatted "peak RAM: X MiB" string. */
    private static String stopCpuPolling() {
        if (cpuPollThread == null) return null;
        cpuPolling = false;
        try { cpuPollThread.join(500); } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
        cpuPollThread = null;

        long peak = peakCpuRssMiB.get();
        if (peak <= 0) return null;
        return "peak RAM: " + peak + " MiB";
    }

    // ── Public API ────────────────────────────────────────────────────────────

    /**
     * Print phase-start header and return a start timestamp.
     * @param label human-readable phase label
     * @param gpu   true if this phase executes on the GPU (reserved for future use)
     * @return System.nanoTime() snapshot for passing to {@link #end}
     */
    public static long begin(String label, boolean gpu) {
        startCpuPolling();
        String tag = gpu ? "[GPU]" : "[CPU]";
        System.err.println();
        System.err.println("  >> " + label + "  " + tag);
        return System.nanoTime();
    }

    /**
     * Print phase-completion line with elapsed time and peak CPU RAM.
     * @param label same label passed to {@link #begin}
     * @param t0    timestamp returned by {@link #begin}
     * @param gpu   true if this phase ran on the GPU
     */
    public static void end(String label, long t0, boolean gpu) {
        long ms = (System.nanoTime() - t0) / 1_000_000;
        String cpuStr  = stopCpuPolling();

        StringBuilder sb = new StringBuilder();
        sb.append("     OK  ").append(ms).append(" ms");
        if (cpuStr != null) {
            sb.append("  |  ").append(cpuStr);
        }
        System.err.println(sb);
    }
}
