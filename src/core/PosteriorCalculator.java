package core;

/**
 * Computes branch length and posterior support 
 */
public class PosteriorCalculator {
    
    private double f1;  // Main topology frequency (supporting species tree)
    private double f2;  // Alternative topology 1 frequency  
    private double f3;  // Alternative topology 2 frequency
    private double n;   // Total effective number of quartets
    private double lambda; // Prior parameter for Beta-Binomial model
    
    private double posterior = -1;
    private double pValue = -1;
    private double branchLength = -1;
    
    private static final double LOG2 = Math.log(2.0);
    private static final double MACHINE_EPSILON = 2.220446049250313e-16;
    private static final boolean DEBUG = false;
    
    /**
     * Constructor for posterior probability calculation
     * 
     * @param f1 Frequency of quartets supporting the main topology
     * @param f2 Frequency of quartets supporting alternative topology 1
     * @param f3 Frequency of quartets supporting alternative topology 2  
     * @param n Total effective number of quartets
     * @param lambda Prior parameter (typically 0.5)
     */
    public PosteriorCalculator(double f1, double f2, double f3, double n, double lambda) {
        this.f1 = f1;
        this.f2 = f2; 
        this.f3 = f3;
        this.n = n;
        this.lambda = lambda;
    }
    
    /**
     * Get posterior probability for the main topology
     */
    public double getPosteriorProbability() {
        if (this.posterior == -1) {
            this.posterior = calculatePosterior();
        }
        return posterior;
    }
    
    /**
     * Get p-value for polytomy test (chi-square test)
     */
    public double getPValue() {
        if (pValue == -1) {
            this.pValue = calculatePValue();
        }
        return pValue;
    }
    
    /**
     * Get branch length based on coalescent theory
     */
    public double getBranchLength() {
        if (branchLength == -1) {
            this.branchLength = calculateBranchLength();
        }
        return branchLength;
    }
    
    /**
     * Calculate posterior probability using Beta-Binomial model
     */
    private double calculatePosterior() {
        if (DEBUG) {
            System.out.println("f1=" + f1 + " f2=" + f2 + " f3=" + f3 + " n=" + n);
            System.out.println("G1: " + G(f1, n));
            System.out.println("G2: " + G(f2, n)); 
            System.out.println("G3: " + G(f3, n));
            System.out.println("r2: " + r(f2));
            System.out.println("r3: " + r(f3));
        }
        
        double g2 = rG(f2);
        double g3 = rG(f3);
        double g1 = G(f1, n);
        
        double denominator = g1 + g2 + g3;
        if (denominator == 0) {
            return 0.0;
        }
        
        double posterior = g1 / denominator;
        
        if (Double.isInfinite(posterior)) {
            throw new RuntimeException("Posterior calculation resulted in infinity. " +
                "Parameters: f1=" + f1 + " f2=" + f2 + " f3=" + f3 + " lambda=" + lambda + " n=" + n);
        }
        
        if (Double.isNaN(posterior)) {
            if (f1 * 3 < n) {
                return 0.0;
            } else {
                throw new RuntimeException("Posterior calculation resulted in NaN. " +
                    "Parameters: f1=" + f1 + " f2=" + f2 + " f3=" + f3 + " lambda=" + lambda + " n=" + n);
            }
        }
        
        return posterior;
    }
    
    /**
     * Calculate p-value using chi-square test for polytomy testing
     */
    private double calculatePValue() {
        if (n <= 15) {
            System.err.println("Warning: Not enough genes (n=" + n + ") for reliable p-value calculation");
            return -2.0; // Indicates insufficient data
        }
        
        double expectedFreq = n / 3.0;
        double chiSquare = Math.pow(f1 - expectedFreq, 2) / expectedFreq +
                          Math.pow(f2 - expectedFreq, 2) / expectedFreq +
                          Math.pow(f3 - expectedFreq, 2) / expectedFreq;
        
        // Chi-square complementary CDF with 2 degrees of freedom
        return chiSquareComplementary(2, chiSquare);
    }
    
    /**
     * Calculate branch length using coalescent theory
     */
    private double calculateBranchLength() {
        return calculateBranchLength(f1, n, lambda);
    }
    
    /**
     * Static method for branch length calculation
     */
    public static double calculateBranchLength(double f, double total, double lambda) {
        double ratio = f / (total + 2 * lambda);
        if (ratio > 1.0/3.0) {
            return -Math.log(1.5 * (1.0 - ratio));
        } else {
            return 0.0;
        }
    }
    
    /**
     * Calculate G function: 1 - incompleteBeta(x+1, n-x+2*lambda, 1/3)
     */
    private double G(double x, double n) {
        double result = 1.0 - incompleteBeta(x + 1, n - x + 2 * lambda, 1.0/3.0);
        if (result <= MACHINE_EPSILON) {
            result = 0.0;
        }
        return result;
    }
    
    /**
     * Calculate r function using beta ratio
     */
    private double r(double mi) {
        double logRatio = logGamma(mi + 1) + logGamma(n - mi + 2 * lambda) -
                         logGamma(f1 + 1) - logGamma(n - f1 + 2 * lambda);
        return Math.exp(LOG2 * (mi - f1) + logRatio);
    }
    
    /**
     * Calculate rG function: G(mi, n) * r(mi) with overflow handling
     */
    private double rG(double mi) {
        double gValue = G(mi, n);
        double rValue = r(mi);
        double result = gValue * rValue;
        
        if (Double.isNaN(result) || Double.isInfinite(result)) {
            if (mi * 3 < n) {
                return 0.0;
            } else if (mi * 3 >= n) {
                if (f1 > mi) {
                    return 0.0;
                } else if (f1 < mi) {
                    return Double.POSITIVE_INFINITY;
                } else {
                    throw new RuntimeException("rG calculation error. " +
                        "Parameters: f1=" + f1 + " f2=" + f2 + " f3=" + f3 + " lambda=" + lambda + " n=" + n);
                }
            } else {
                return Double.POSITIVE_INFINITY;
            }
        }
        
        return result;
    }
    
    /**
     * Log-gamma function approximation using Stirling's approximation for large values
     * and series expansion for smaller values
     */
    private double logGamma(double x) {
        if (x <= 0) {
            throw new IllegalArgumentException("logGamma: x must be positive");
        }
        
        // For large x, use Stirling's approximation
        if (x > 12) {
            double logSqrt2Pi = 0.9189385332046727; // log(sqrt(2*pi))
            return logSqrt2Pi + (x - 0.5) * Math.log(x) - x + 1.0/(12.0*x) - 1.0/(360.0*x*x*x);
        }
        
        // For smaller x, use recurrence relation and series
        double result = 0.0;
        while (x < 8) {
            result -= Math.log(x);
            x += 1.0;
        }
        
        // Now use Stirling's approximation for x >= 8
        double logSqrt2Pi = 0.9189385332046727;
        result += logSqrt2Pi + (x - 0.5) * Math.log(x) - x + 1.0/(12.0*x);
        
        return result;
    }
    
    /**
     * Incomplete beta function I_x(a,b) using continued fraction
     */
    private double incompleteBeta(double a, double b, double x) {
        if (x < 0 || x > 1) {
            throw new IllegalArgumentException("x must be between 0 and 1");
        }
        if (x == 0 || x == 1) {
            return x;
        }
        
        // Use symmetry relation if needed for better convergence
        boolean symmetric = false;
        if (x > (a + 1.0) / (a + b + 2.0)) {
            symmetric = true;
            double temp = a;
            a = b;
            b = temp;
            x = 1.0 - x;
        }
        
        double result = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + 
                               a * Math.log(x) + b * Math.log(1.0 - x)) / a;
        
        // Continued fraction evaluation
        double cf = continuedFractionBeta(a, b, x);
        result *= cf;
        
        return symmetric ? 1.0 - result : result;
    }
    
    /**
     * Continued fraction for incomplete beta function
     */
    private double continuedFractionBeta(double a, double b, double x) {
        final int maxIterations = 200;
        final double epsilon = 1e-15;
        
        double am = 1.0;
        double bm = 1.0;
        double az = 1.0;
        double qab = a + b;
        double qap = a + 1.0;
        double qam = a - 1.0;
        double bz = 1.0 - qab * x / qap;
        
        for (int m = 1; m <= maxIterations; m++) {
            double em = m;
            double tem = em + em;
            double d = em * (b - m) * x / ((qam + tem) * (a + tem));
            double ap = az + d * am;
            double bp = bz + d * bm;
            d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));
            double app = ap + d * az;
            double bpp = bp + d * bz;
            double aold = az;
            am = ap / bpp;
            bm = bp / bpp;
            az = app / bpp;
            bz = 1.0;
            
            if (Math.abs(az - aold) < epsilon * Math.abs(az)) {
                return az;
            }
        }
        
        throw new RuntimeException("Continued fraction did not converge in incompleteBeta");
    }
    
    /**
     * Chi-square complementary cumulative distribution function
     * P(X > x) where X ~ Chi-square(df)
     */
    private double chiSquareComplementary(int degreesOfFreedom, double x) {
        if (x <= 0) return 1.0;
        if (degreesOfFreedom <= 0) throw new IllegalArgumentException("Degrees of freedom must be positive");
        
        // For df=2, we have exponential distribution: P(X > x) = exp(-x/2)
        if (degreesOfFreedom == 2) {
            return Math.exp(-x / 2.0);
        }
        
        // For general case, use incomplete gamma function
        // P(X > x) = 1 - P(X <= x) = 1 - (1/Gamma(df/2)) * gamma(df/2, x/2)
        // where gamma(a,x) is the lower incomplete gamma function
        double a = degreesOfFreedom / 2.0;
        return 1.0 - incompleteGamma(a, x / 2.0);
    }
    
    /**
     * Incomplete gamma function gamma(a,x)/Gamma(a)
     */
    private double incompleteGamma(double a, double x) {
        if (x < 0) return 0.0;
        if (x == 0) return 0.0;
        
        if (x < a + 1.0) {
            // Use series expansion
            return incompleteGammaSeries(a, x);
        } else {
            // Use continued fraction
            return 1.0 - incompleteGammaContinuedFraction(a, x);
        }
    }
    
    /**
     * Incomplete gamma function using series expansion
     */
    private double incompleteGammaSeries(double a, double x) {
        final int maxIterations = 200;
        final double epsilon = 1e-15;
        
        double sum = 1.0 / a;
        double term = 1.0 / a;
        
        for (int n = 1; n <= maxIterations; n++) {
            term *= x / (a + n);
            sum += term;
            if (Math.abs(term) < epsilon * Math.abs(sum)) {
                break;
            }
        }
        
        return Math.exp(-x + a * Math.log(x) - logGamma(a)) * sum;
    }
    
    /**
     * Incomplete gamma function using continued fraction
     */
    private double incompleteGammaContinuedFraction(double a, double x) {
        final int maxIterations = 200;
        final double epsilon = 1e-15;
        
        double b = x + 1.0 - a;
        double c = 1e30;
        double d = 1.0 / b;
        double h = d;
        
        for (int i = 1; i <= maxIterations; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (Math.abs(d) < epsilon) d = epsilon;
            c = b + an / c;
            if (Math.abs(c) < epsilon) c = epsilon;
            d = 1.0 / d;
            double del = d * c;
            h *= del;
            if (Math.abs(del - 1.0) < epsilon) {
                break;
            }
        }
        
        return Math.exp(-x + a * Math.log(x) - logGamma(a)) * h;
    }
    
    @Override
    public String toString() {
        return String.format("Posterior: %.6f", getPosteriorProbability());
    }
    
    /**
     * Detailed string representation
     */
    public String toDetailedString() {
        return String.format("f1=%.2f, f2=%.2f, f3=%.2f, n=%.2f, posterior=%.6f, branchLength=%.6f", 
                           f1, f2, f3, n, getPosteriorProbability(), getBranchLength());
    }
}
