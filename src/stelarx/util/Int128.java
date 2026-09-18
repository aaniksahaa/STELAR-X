package stelarx.util;

import java.math.BigInteger;

/**
 * Minimal immutable signed 128-bit integer, represented as a high signed 64-bit
 * word and a low <em>unsigned</em> 64-bit word ({@code value = hi · 2^64 + (lo & 0xFFFF…)}).
 *
 * <p>Used for exact quartet-score accumulation when the value exceeds the range
 * of {@code long} (very large taxon sets).  Only the operations the scoring +
 * inference DP actually need are provided, each implemented with a handful of
 * primitive instructions (no {@code BigInteger} on the hot path):
 * {@link #add}, {@link #compareTo}, {@link #halve}, and the small-operand
 * builders {@link #ofLong} / {@link #mulLong}.
 *
 * <p>All scores in STELAR-X are non-negative; the signed representation is used
 * only so the DP can start from a {@code null}/lowest sentinel and compare.
 * {@link #toString} (rare, for logging) goes through {@code BigInteger}.
 */
public final class Int128 {

    /** Constant zero. */
    public static final Int128 ZERO = new Int128(0L, 0L);

    public final long hi;   // signed high word
    public final long lo;   // unsigned low word

    public Int128(long hi, long lo) {
        this.hi = hi;
        this.lo = lo;
    }

    /** Sign-extend a long into 128 bits. */
    public static Int128 ofLong(long v) {
        return new Int128(v < 0 ? -1L : 0L, v);
    }

    /**
     * Exact 64×64→128 product of two <em>non-negative</em> longs.
     * Uses {@link Math#multiplyHigh} for the upper 64 bits (single instruction
     * on modern JITs); the low word is the wrapping product.
     */
    public static Int128 mulLong(long a, long b) {
        long lo = a * b;
        long hi = Math.multiplyHigh(a, b);
        return new Int128(hi, lo);
    }

    /** this + other (128-bit add with carry). */
    public Int128 add(Int128 other) {
        long rlo = this.lo + other.lo;
        // carry out of the low word when the unsigned sum wrapped below an addend
        long carry = (Long.compareUnsigned(rlo, this.lo) < 0) ? 1L : 0L;
        long rhi = this.hi + other.hi + carry;
        return new Int128(rhi, rlo);
    }

    /**
     * Multiply this non-negative value by a small non-negative long scalar.
     * Correct as long as the true product fits in 128 bits (guaranteed for the
     * score magnitudes here: |this| ≤ ~2^72, scalar ≤ ~2^14).
     */
    public Int128 mulScalar(long f) {
        // lo is an unsigned limb and may therefore be negative as a Java long.
        // Math.multiplyHigh/mulLong would sign-extend it and lose f*2^64 whenever
        // bit 63 is set. Use the unsigned high half for this limb product.
        long productLo = this.lo * f;
        long productHi = Math.unsignedMultiplyHigh(this.lo, f);
        long hiAdd = this.hi * f;                // high-word contribution (fits here)
        return new Int128(productHi + hiAdd, productLo);
    }

    /** Logical right shift by one (exact divide-by-two for non-negative values). */
    public Int128 halve() {
        long newLo = (lo >>> 1) | (hi << 63);
        long newHi = hi >>> 1;
        return new Int128(newHi, newLo);
    }

    /** Signed 128-bit comparison. */
    public int compareTo(Int128 o) {
        if (this.hi != o.hi) return Long.compare(this.hi, o.hi);
        return Long.compareUnsigned(this.lo, o.lo);
    }

    /** Approximate floating-point value (for stats/double consumers). */
    public double toDouble() {
        // hi·2^64 + lo(unsigned).  Sufficient precision for reporting/conversion.
        double loD = (lo >>> 1) * 2.0 + (lo & 1L);   // unsigned long → double
        return (double) hi * 18446744073709551616.0 /* 2^64 */ + loD;
    }

    /** Exact value as BigInteger (used by {@link #toString}; not hot-path). */
    public BigInteger toBigInteger() {
        BigInteger hiB = BigInteger.valueOf(hi).shiftLeft(64);
        BigInteger loB = new BigInteger(Long.toUnsignedString(lo));
        return hiB.add(loB);
    }

    @Override public String toString() { return toBigInteger().toString(); }

    @Override public boolean equals(Object o) {
        if (!(o instanceof Int128)) return false;
        Int128 x = (Int128) o;
        return hi == x.hi && lo == x.lo;
    }

    @Override public int hashCode() { return Long.hashCode(hi) * 31 + Long.hashCode(lo); }
}
