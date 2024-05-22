const std = @import("std");
const testing = std.testing;
const builtin = @import("builtin");
// A fast generator of discrete, bounded
// [Zipf-distributed](https://en.wikipedia.org/wiki/Zipf's_law) random numbers.
//
// For a random variable `X` whose values are distributed according to this distribution, the
// probability mass function is given by
//
// ```ignore
// P(X = k) = H(N,s) * 1 / k^s for k = 1,2,...,N
// ```
//
// `H(N,s)` is the normalizing constant which corresponds to the generalized harmonic number
// of order `N` of `s`.
//
// # Example
//
// ```
// use rand::distributions::Distribution;
//
// let mut rng = rand::thread_rng();
// let mut zipf = zipf::ZipfDistribution::new(1000, 1.03).unwrap();
// let sample = zipf.sample(&mut rng);
// ```
//
// This implementation is effectively a direct port of Apache Common's
// [RejectionInversionZipfSampler][ref], written in Java. It is based on the method described by
// Wolfgang Hörmann and Gerhard Derflinger in [*Rejection-inversion to generate variates from
// monotone discrete distributions*](https://dl.acm.org/citation.cfm?id=235029) from *ACM
// Transactions on Modeling and Computer Simulation (TOMACS) 6.3 (1996)*.
//
// [ref]: https://github.com/apache/commons-rng/blob/6a1b0c16090912e8fc5de2c1fb5bd8490ac14699/commons-rng-sampling/src/main/java/org/apache/commons/rng/sampling/distribution/RejectionInversionZipfSampler.java
//

// Random number generator that generates Zipf-distributed random numbers using rejection

/// inversion.
///
const ZipfDistributionError = error{
    ElementsZero,
    ExponentZeroOrNegative,
};
pub const ZipfDistribution = struct {
    // number elemnts
    num_elements: f64,
    // exponent parameter of the distribution
    exponent: f64,
    // hIntegral(1.5) - 1)
    h_integral_x1: f64,
    // hIntegral(num_elements + 0.5)
    h_integral_num_elements: f64,
    // `2 - hIntegralInverse(hIntegral(2.5) - h(2)}`
    s: f64,

    // Constructor
    pub fn init(num_elements: u64, exponent: f64) !ZipfDistribution {
        if (num_elements == 0) {
            return ZipfDistributionError.ElementsZero;
        }
        if (exponent <= 0.0) {
            return ZipfDistributionError.ExponentZeroOrNegative;
        }
        const num_elements_float: f64 = @floatFromInt(num_elements);

        return ZipfDistribution{
            // used for calculation later therefore internally as float
            .num_elements = num_elements_float,
            .exponent = exponent,
            .h_integral_x1 = ZipfDistribution.h_integral(1.5, exponent) - 1.0,
            .h_integral_num_elements = ZipfDistribution.h_integral(num_elements_float + 0.5, exponent),
            .s = 2.0 - ZipfDistribution.h_integral_inv(
                ZipfDistribution.h_integral(2.5, exponent) - ZipfDistribution.h(2.0, exponent),
                exponent,
            ),
        };
    }

    pub fn next(self: *ZipfDistribution, rand: *const std.rand.Random) u64 {
        const hnum = self.h_integral_num_elements;

        while (true) {
            const u: f64 = hnum + rand.float(f64) * (self.h_integral_x1 - hnum);
            const x: f64 = ZipfDistribution.h_integral_inv(u, self.exponent);

            // Limit k to the range [1, num_elements] if it would be outside
            // due to numerical inaccuracies.
            const k64 = @min(@max(x, 1.0), self.num_elements);
            // float -> integer rounds towards zero, so we add 0.5
            // to prevent bias towards k == 1
            const k = @max(1, (k64 + 0.5));

            // Here, the distribution of k is given by:
            //
            //   P(k = 1) = C * (hIntegral(1.5) - h_integral_x1) = C
            //   P(k = m) = C * (hIntegral(m + 1/2) - hIntegral(m - 1/2)) for m >= 2
            //
            // where C = 1 / (h_integral_num_elements - h_integral_x1)
            if (k64 - x <= self.s or u >= (ZipfDistribution.h_integral(k64 + 0.5, self.exponent)) - ZipfDistribution.h(k64, self.exponent)) {

                // Case k = 1:
                //
                //   The right inequality is always true, because replacing k by 1 gives
                //   u >= hIntegral(1.5) - h(1) = h_integral_x1 and u is taken from
                //   (h_integral_x1, h_integral_num_elements].
                //
                //   Therefore, the acceptance rate for k = 1 is P(accepted | k = 1) = 1
                //   and the probability that 1 is returned as random value is
                //   P(k = 1 and accepted) = P(accepted | k = 1) * P(k = 1) = C = C / 1^exponent
                //
                // Case k >= 2:
                //
                //   The left inequality (k - x <= s) is just a short cut
                //   to avoid the more expensive evaluation of the right inequality
                //   (u >= hIntegral(k + 0.5) - h(k)) in many cases.
                //
                //   If the left inequality is true, the right inequality is also true:
                //     Theorem 2 in the paper is valid for all positive exponents, because
                //     the requirements h'(x) = -exponent/x^(exponent + 1) < 0 and
                //     (-1/hInverse'(x))'' = (1+1/exponent) * x^(1/exponent-1) >= 0
                //     are both fulfilled.
                //     Therefore, f(x) = x - hIntegralInverse(hIntegral(x + 0.5) - h(x))
                //     is a non-decreasing function. If k - x <= s holds,
                //     k - x <= s + f(k) - f(2) is obviously also true which is equivalent to
                //     -x <= -hIntegralInverse(hIntegral(k + 0.5) - h(k)),
                //     -hIntegralInverse(u) <= -hIntegralInverse(hIntegral(k + 0.5) - h(k)),
                //     and finally u >= hIntegral(k + 0.5) - h(k).
                //
                //   Hence, the right inequality determines the acceptance rate:
                //   P(accepted | k = m) = h(m) / (hIntegrated(m+1/2) - hIntegrated(m-1/2))
                //   The probability that m is returned is given by
                //   P(k = m and accepted) = P(accepted | k = m) * P(k = m)
                //                         = C * h(m) = C / m^exponent.
                //
                // In both cases the probabilities are proportional to the probability mass
                // function of the Zipf distribution.
                return @intFromFloat(k);
            }
        }
    }

    /// Computes `H(x)`, defined as
    ///
    ///  - `(x^(1 - exponent) - 1) / (1 - exponent)`, if `exponent != 1`
    ///  - `log(x)`, if `exponent == 1`
    ///
    /// `H(x)` is an integral function of `h(x)`, the derivative of `H(x)` is `h(x)`.
    fn h_integral(x: f64, exponent: f64) f64 {
        const log_x = @log(x);
        return helper2((1.0 - exponent) * log_x) * log_x;
    }

    /// Computes `h(x) = 1 / x^exponent`
    fn h(x: f64, exponent: f64) f64 {
        return std.math.exp(-exponent * @log(x));
    }

    /// The inverse function of `H(x)`.
    /// Returns the `y` for which `H(y) = x`.
    fn h_integral_inv(x: f64, exponent: f64) f64 {
        var t: f64 = x * (1.0 - exponent);
        if (t < -1.0) {
            // Limit value to the range [-1, +inf).
            // t could be smaller than -1 in some rare cases due to numerical errors.
            t = -1.0;
        }
        return std.math.exp((helper1(t) * x));
    }
    /// Helper function to calculate `(exp(x) - 1) / x`.
    /// A Taylor series expansion is used, if x is close to 0.
    fn helper2(x: f64) f64 {
        if (@abs(x) > 1e-8) {
            return std.math.expm1(x) / x;
        } else {
            return 1.0 + x * 0.5 * (1.0 + x * 1.0 / 3.0 * (1.0 + 0.25 * x));
        }
    }

    /// Helper function that calculates `log(1 + x) / x`.
    /// A Taylor series expansion is used, if x is close to 0.
    fn helper1(x: f64) f64 {
        if (@abs(x) > 1e-8) {
            return std.math.log1p(x) / x;
        } else {
            return 1.0 - x * (0.5 - x * (1.0 / 3.0 - 0.25 * x));
        }
    }
};

const expect = std.testing.expect;
test {
    // To run nested container tests, either, call `refAllDecls` which will
    // reference all declarations located in the given argument.
    // `@This()` is a builtin function that returns the innermost container it is called from.
    // In this example, the innermost container is this file (implicitly a struct).
    std.testing.refAllDecls(@This());

    // or, reference each container individually from a top-level test declaration.
    // The `_ = C;` syntax is a no-op reference to the identifier `C`.
    _ = S;
}

const S = struct {
    test "one" {
        try test_helper(1.00);
    }

    test "two" {
        try test_helper(2.00);
    }

    test "three" {
        try test_helper(3.00);
    }
    test "float" {
        try test_helper(1.08);
    }

    test "errors" {
        {
            const zipf = ZipfDistribution.init(0, 1.0);
            try std.testing.expectError(ZipfDistributionError.ElementsZero, zipf);
        }
        {
            const zipf = ZipfDistribution.init(100, 0.0);
            try std.testing.expectError(ZipfDistributionError.ExponentZeroOrNegative, zipf);
        }
        {
            const zipf = ZipfDistribution.init(100, -1.0);
            try std.testing.expectError(ZipfDistributionError.ExponentZeroOrNegative, zipf);
        }
    }

    fn test_helper(alpha: f64) !void {
        const N: usize = 100;

        // as the alpha increases, we need more samples, since the frequency in the tail grows so low

        const samples: usize = @intFromFloat((std.math.pow(f64, 2.0, alpha) * 5000000.0));

        var prng = std.rand.DefaultPrng.init(blk: {
            var seed: u64 = undefined;
            try std.posix.getrandom(std.mem.asBytes(&seed));
            break :blk seed;
        });
        const rand = prng.random();

        var zipf_distribution = try ZipfDistribution.init(N, alpha);
        var harmonic: f64 = 0.0;
        for (1..N) |n| {
            const n_f: f64 = @floatFromInt(n);
            harmonic += 1.0 / std.math.pow(f64, n_f, alpha);
        }

        // sample a bunch
        var buckets: [N]usize = std.mem.zeroes([N]usize);
        for (0..samples) |_| {
            const sample = zipf_distribution.next(&rand);
            buckets[sample - 1] += 1;
        }

        // for each bucket, see that frequency roughly matches what we expect
        // note that the ratios here are ratios _of fractions_, so 0.5 does not mean we're off by
        // 50%, it means we're off by 50% _of the frequency_. in other words, if the frequency was
        // supposed to be 0.10, and it was actually 0.05, that's a 50% difference.
        for (1..N) |i| {
            const b_i: f64 = @floatFromInt(buckets[i]);
            const freq: f64 = b_i / @as(f64, @floatFromInt(samples));
            const expected = (1.0 / std.math.pow(f64, @as(f64, @floatFromInt(i)) + 1.0, alpha)) / harmonic;

            const off_by = @abs(expected - freq);
            try expect(off_by < 0.1); // never be off by more than 10% in absolute terms

            const good: bool = (off_by < expected);
            if (!good) {
                std.log.err("got {}, expected {} for k = {}", .{ freq, expected, i + 1 });
            }
        }
    }
};
