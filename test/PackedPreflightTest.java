import stelarx.completion.SimilarityMatrixBuilder;

import java.util.List;

/** Ensures an undersized heap fails before any multi-GiB matrix allocation. */
public final class PackedPreflightTest {
    public static void main(String[] args) {
        try {
            SimilarityMatrixBuilder.buildCPU(List.of(), 50_000);
            throw new AssertionError("50k packed preflight unexpectedly accepted a 1 GiB heap");
        } catch (IllegalArgumentException expected) {
            if (!expected.getMessage().contains("Exact packed similarity accumulators")
                    || !expected.getMessage().contains("Increase --xmx")) {
                throw expected;
            }
        }
        System.out.println("Packed matrix heap preflight: PASS");
    }
}
