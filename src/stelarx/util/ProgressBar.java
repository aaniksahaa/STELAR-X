package stelarx.util;

import stelarx.Logging;
import java.util.concurrent.atomic.AtomicLong;

/**
 * tqdm-style inline progress bar for CPU work loops, written to stderr.
 *
 * Format (TTY):
 *   ▸  Parsing trees  [████████████░░░░░░░░░░░░░░░░]  500/1000 (50%)  [00:02<00:02, 231it/s]
 *
 * Always active at INFO level — prints unconditionally regardless of TTY,
 * matching the behaviour of the native GPU progress output.
 *
 * Thread-safe: update() may be called from multiple threads concurrently;
 * a CAS-based rate-limiter ensures at most one repaint per INTERVAL_NS.
 */
public class ProgressBar {

    private static final boolean COLOR = stelarx.Banner.useColor();

    private static final String RST  = "\033[0m";
    private static final String DIM  = "\033[2m";
    private static final String CYAN = "\033[36m";
    private static final String YLW  = "\033[33m";

    private static final int  BAR_WIDTH   = 28;
    private static final long INTERVAL_NS = 150_000_000L; // 150 ms

    private final String     label;
    private final int        total;
    private final long       startNs;
    private final AtomicLong lastPaintNs = new AtomicLong(0L);

    // -------------------------------------------------------------------------

    public ProgressBar(String label, int total) {
        this.label   = label;
        this.total   = total;
        this.startNs = System.nanoTime();
        if (active() && total > 0) {
            render(0);
            lastPaintNs.set(System.nanoTime());
        }
    }

    /**
     * Report progress. Safe to call from multiple threads.
     * @param done number of items completed so far (1-based)
     */
    public void update(int done) {
        if (!active() || total <= 0) return;
        long now  = System.nanoTime();
        long last = lastPaintNs.get();
        if (now - last < INTERVAL_NS && done < total) return;
        if (!lastPaintNs.compareAndSet(last, now)) return;
        render(done);
    }

    /** Call once after the loop finishes. Prints the completed bar on its own line. */
    public void done() {
        if (!active() || total <= 0) return;
        render(total);
        System.err.println();
    }

    // -------------------------------------------------------------------------

    private boolean active() {
        return Logging.getLevel() >= Logging.INFO;
    }

    private void render(int done) {
        int pct    = (int)(100L * done / total);
        int filled = (int)(1L   * BAR_WIDTH * done / total);

        double elapsedSec = (System.nanoTime() - startNs) / 1e9;
        double rate       = (done > 0 && elapsedSec > 0) ? done / elapsedSec : 0;
        double etaSec     = (rate > 0 && done < total)   ? (total - done) / rate : 0;

        StringBuilder sb = new StringBuilder("     ");
        sb.append(COLOR ? DIM + "▸  " + RST : "▸  ");
        sb.append(label).append("  ");
        // bar
        sb.append(COLOR ? CYAN : "");
        sb.append('[');
        for (int i = 0; i < BAR_WIDTH; i++) sb.append(i < filled ? '█' : '░');
        sb.append(']');
        if (COLOR) sb.append(RST);
        // count + pct
        sb.append(String.format("  %d/%d (%d%%)", done, total, pct));
        // tqdm-style timing block
        sb.append("  ");
        sb.append(COLOR ? DIM : "");
        sb.append('[');
        sb.append(fmtTime(elapsedSec));
        sb.append('<');
        sb.append(done > 0 ? fmtTime(etaSec) : "?");
        sb.append(", ");
        if (COLOR) sb.append(RST + YLW);
        sb.append(fmtRate(rate));
        if (COLOR) sb.append(RST + DIM);
        sb.append(']');
        if (COLOR) sb.append(RST);

        System.err.print(sb + "\r");
    }

    // ── Formatting helpers ────────────────────────────────────────────────────

    private static String fmtTime(double sec) {
        long s = (long) sec;
        long m = s / 60; s %= 60;
        if (m < 60) return String.format("%02d:%02d", m, s);
        long h = m / 60; m %= 60;
        return String.format("%d:%02d:%02d", h, m, s);
    }

    private static String fmtRate(double rate) {
        if (rate <= 0)   return "?it/s";
        if (rate >= 1.0) return String.format("%.1fit/s", rate);
        // very slow: show seconds-per-item instead
        return String.format("%.2fs/it", 1.0 / rate);
    }
}
