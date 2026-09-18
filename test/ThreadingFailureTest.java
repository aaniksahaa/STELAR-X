import stelarx.util.Threading;

/** Proves that parallel worker failures reach the caller instead of being lost. */
public final class ThreadingFailureTest {
    public static void main(String[] args) {
        Threading.start(4);
        try {
            try {
                Threading.processRangeParallel(100, i -> {
                    if (i == 37) throw new IllegalStateException("worker-marker-37");
                });
                throw new AssertionError("parallel worker failure was swallowed");
            } catch (IllegalStateException expected) {
                if (!"worker-marker-37".equals(expected.getMessage())) throw expected;
            }
        } finally {
            Threading.shutdown();
        }
        System.out.println("Threading failure propagation: PASS");
    }
}
