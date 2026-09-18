package stelarx.weight;

/** Boundary/monotonicity tests for LONG versus INT128/DOUBLE dispatch. */
public final class WeightModeBoundaryTest {
    private static void check(boolean condition, String message) {
        if (!condition) throw new AssertionError(message);
    }

    public static void main(String[] args) {
        String forced = args.length == 0 ? "none" : args[0];
        if (forced.equals("double")) {
            check(WeightTable.needsDoubleAccumulation(3, 1),
                "STELARX_WEIGHT_FORCE_DOUBLE did not force the wide path");
            System.out.println("Weight mode forced-double dispatch: PASS");
            return;
        }
        if (forced.equals("long")) {
            check(!WeightTable.needsDoubleAccumulation(Integer.MAX_VALUE, Integer.MAX_VALUE),
                "STELARX_WEIGHT_FORCE_LONG did not force the long path");
            System.out.println("Weight mode forced-long dispatch: PASS");
            return;
        }

        check(!WeightTable.needsDoubleAccumulation(3, 1), "smallest score should fit long");
        check(!WeightTable.needsDoubleAccumulation(1_000, 1_000),
            "ordinary problem unexpectedly selected a wide score type");
        check(WeightTable.needsDoubleAccumulation(2_000_000, 1_000_000),
            "large problem did not select a wide score type");

        for (int genes : new int[] {1, 10, 1_000, 1_000_000}) {
            boolean seenWide = false;
            int transitions = 0;
            for (int n = 3; n <= 4_000_000; n = next(n)) {
                boolean wide = WeightTable.needsDoubleAccumulation(n, genes);
                if (wide && !seenWide) {
                    seenWide = true;
                    transitions++;
                }
                check(!seenWide || wide, "dispatch ceased being monotonic for genes=" + genes);
            }
            check(transitions <= 1, "multiple precision transitions for genes=" + genes);
        }
        System.out.println("Weight numeric-mode boundaries: PASS");
    }

    private static int next(int n) {
        long candidate = Math.max((long)n + 1L, (long)Math.ceil(n * 1.17));
        return candidate > 4_000_000L ? 4_000_001 : (int)candidate;
    }
}
