package stelarx.util;

import java.math.BigInteger;
import java.util.Random;

/** Differential arithmetic tests against BigInteger for the exact score type. */
public final class Int128Test {
    private static final BigInteger MASK64 = BigInteger.ONE.shiftLeft(64).subtract(BigInteger.ONE);

    private static Int128 fromBigInteger(BigInteger value) {
        return new Int128(value.shiftRight(64).longValue(), value.and(MASK64).longValue());
    }

    private static void equal(BigInteger expected, Int128 actual, String operation) {
        if (!expected.equals(actual.toBigInteger())) {
            throw new AssertionError(operation + ": expected=" + expected
                + " actual=" + actual.toBigInteger());
        }
    }

    public static void main(String[] args) {
        long[] edges = {0L, 1L, 2L, 3L, Integer.MAX_VALUE, 1L << 32,
                        Long.MAX_VALUE - 1, Long.MAX_VALUE};
        for (long value : edges) {
            equal(BigInteger.valueOf(value), Int128.ofLong(value), "ofLong");
            equal(BigInteger.valueOf(value).shiftRight(1),
                Int128.ofLong(value).halve(), "halve-edge");
        }

        Random random = new Random(0x1285E1A7L);
        int checks = edges.length * 2;
        for (int i = 0; i < 20_000; i++) {
            long a = random.nextLong() & Long.MAX_VALUE;
            long b = random.nextLong() & Long.MAX_VALUE;
            BigInteger product = BigInteger.valueOf(a).multiply(BigInteger.valueOf(b));
            equal(product, Int128.mulLong(a, b), "mulLong");

            BigInteger x = new BigInteger(120, random);
            BigInteger y = new BigInteger(120, random);
            equal(x.add(y), fromBigInteger(x).add(fromBigInteger(y)), "add");
            equal(x.shiftRight(1), fromBigInteger(x).halve(), "halve");

            BigInteger smallBase = new BigInteger(88, random);
            long scalar = random.nextInt(1_000_001);
            equal(smallBase.multiply(BigInteger.valueOf(scalar)),
                fromBigInteger(smallBase).mulScalar(scalar), "mulScalar");

            int expectedComparison = x.compareTo(y);
            int actualComparison = fromBigInteger(x).compareTo(fromBigInteger(y));
            if (Integer.signum(expectedComparison) != Integer.signum(actualComparison)) {
                throw new AssertionError("compareTo mismatch: " + x + " vs " + y);
            }
            checks += 5;
        }
        System.out.println("Int128/BigInteger differential arithmetic: PASS (" + checks
            + " checks)");
    }
}
