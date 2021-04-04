// CORDIC Demonstration in Rust
//
// Feel free to copy-paste into https://play.rust-lang.org
//
// This uses fixed-point integers to demonstrate the efficiency
// realized by power-of-two multiplications and divisions.
//
// Floating point integers are normally broken up into three segments
//   1. Base
//   2. Exponent
//   3. "Special" bits
//
// Floating point works similarly to scientific notation, where base is
// the base, exponent is the exponent (although normally expressed as
// 2^n instead of 10^n). The "Special" bits define things like
// Not-a-Number (NaN) for cases like 0/0, +Inf and -Inf in cases
// where the base and exponents do not have enough bits to express
// the full value
//
// Fixed point is just an integer with a set number of places reserved
// at the end for fractional components
//
// If you have any questions, feel free to email me at djh4@illinois.edu
//
// NOTE: This entire file is 100% self contained and will work on
// play.rust-lang.org (you may need to replace theta and iters with
// constants, but that should be it)

use std::cmp::{Ordering, PartialEq, PartialOrd};
use std::fmt;
use std::ops::{Add, Div, Mul, Rem, Sub};

// FixedPoint actually wraps floating point numbers, so currently there
// isn't a difference, but this is an opportunity to come up with your own
// real number representation. External libraries do exist that can
// handle this
pub struct FixedPoint {
    val: f64,
}

impl FixedPoint {
    fn new(val: f64) -> Self {
        Self { val }
    }
}

impl Copy for FixedPoint {}

impl Clone for FixedPoint {
    fn clone(&self) -> Self {
        Self { val: self.val }
    }
}

impl Add for FixedPoint {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            val: self.val + other.val,
        }
    }
}

impl Sub for FixedPoint {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            val: self.val - other.val,
        }
    }
}

impl Mul for FixedPoint {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            val: self.val * other.val,
        }
    }
}

impl Div for FixedPoint {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self {
            val: self.val / other.val,
        }
    }
}

impl Rem for FixedPoint {
    type Output = Self;
    fn rem(self, modulus: FixedPoint) -> Self {
        Self {
            val: self.val.rem(modulus.val),
        }
    }
}

impl PartialOrd for FixedPoint {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.val.partial_cmp(&other.val)
    }
}

impl PartialEq for FixedPoint {
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}

impl fmt::Display for FixedPoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}

fn cordic(mut theta: FixedPoint, iters: usize) -> [FixedPoint; 2] {
    // CORDIC *can* calculate theta even if theta wraps around, but
    // at the expense of accuracy. We make sure that theta is in
    // its lowest valid representation (remainder after division
    // by 2*pi). Another issue with floating point numbers is
    // buildup of error across iterations. Fixed point arithmetic
    // can fix this, but I haven't gotten around to implementing that.
    //
    // NOTE: The "buildup of error across iterations" is why the
    // current unit test fails. Moving over to proper fixed points
    // should fix this
    theta = theta.rem(FixedPoint::new(2.0 * std::f64::consts::PI));
    // CORDIC (for trig functions, at least) does require some
    // compile time constants. However, this is far more space
    // efficient than naively storing sine itself. These are
    // included here to allow for larger iteration countes, but
    // any reasonable implementation would pre-compute a certain
    // amount and keep them in a global array
    // atan(2^-x)
    let angles = (0..(iters as i32 + 2))
        .map(|x| FixedPoint::new((2_f64.powi(-1 * x)).atan()))
        .collect::<Vec<FixedPoint>>();

    // cumprod(1 / sqrt(1 + 2^-2y))
    // NOTE The cumulative product is done by re-calculating and multiplying
    // all elements of the vector together with fold()
    //
    // NOTE 2: This is static for a given number of iters, so in instances
    // where we only compute a set number of iterations, this can be
    // computed ahead of time
    let kvalue = FixedPoint::new(
        (0_i32..(iters as i32 - 1))
            .map(|y| 1.0_f64 / (1.0_f64 + 2_f64.powi(-2 * y)).sqrt().abs())
            .fold(1.0_f64, |x, y| x * y),
    );

    let mut poweroftwo = FixedPoint::new(1.0);
    let mut angle = angles[0];
    let fixed_point_zero = FixedPoint::new(0.0);
    let fixed_point_pos_one = FixedPoint::new(1.0);
    let fixed_point_neg_one = FixedPoint::new(-1.0);
    let fixed_point_pos_two = FixedPoint::new(2.0);

    let mut v = [fixed_point_pos_one, fixed_point_zero]; // Initialize as cos = 1, sine = 0
    for i in 0..iters {
        let sigma_is_neg = theta < fixed_point_zero;

        // v = R * v
        // NOTE: Matrix is always of the form
        // [ 1.0, -factor; factor, 1.0 ]
        //
        // You can imagine `v` as a vector whose cosine is
        // one and sine is zero. The following matrix is a
        // rotation matrix, and it is normally of the form
        // [ cos theta, -sin theta; sin theta, cos theta ].
        // However, this also magnfiies the vector by some
        // set amount, and we account for all these magnitude
        // increases in one multiplication at the end of the
        // calculation
        //
        // The following simplifies down to a rotation of
        // tan^-1(2^-i) and a increase in magnitude of
        // (1 + 2^(-2j))^(1/2)

        let factor = if sigma_is_neg {
            // NOTE: Almost all compilers will optimize multiplication by
            // -1 to be flipping a single bit, so for performance reasons
            // this isn't a normal multiplication
            //
            // HOWEVER, if -1 were stored in a variable, that optimization
            // isn't as obvious, so the compiler might not catch it. We store
            // sigma_is_neg instead of sigma = -1 or sigma = 1 because of this
            FixedPoint::new(-1.0) * poweroftwo
        } else {
            poweroftwo
        };
        let matrix = [
            [fixed_point_pos_one, fixed_point_neg_one * factor],
            [factor, fixed_point_pos_one],
        ];

        // NOTE: Every element of the matrix is either 1 or 2^-n for some
        // number n. An efficient implementation would directly apply
        // those bit-shifts here, but for the sake of explaining this,
        // we will multiply the two values without any optimizations
        v = [
            v[0] + fixed_point_neg_one * factor * v[1],
            factor * v[0] + v[1],
        ];

        // sigma
        theta = if sigma_is_neg {
            theta + angle
        } else {
            theta - angle
        };
        poweroftwo = poweroftwo / fixed_point_pos_two;

        angle = angles[i as usize + 1];
    }

    // Scale vector back such that magnitude is 1
    // NOTE: This can be done either by keeping track of the
    // initial values or performing a square root. If the machine
    // is slow enough that CORDIC is practical (i.e. expensive
    // hardware multiplication), then it is too slow for square roots
    // and divisions
    [v[0] * kvalue, v[1] * kvalue]
}

fn taylor(mut theta: FixedPoint, iters: usize) -> [FixedPoint; 2] {
    todo!()
}


fn main() {
    // Pull parameters from string, should be called as either
    // ./cordic-rs [theta] [iters]
    // or
    // cargo run [theta] [iters]
    let mode = std::env::args().nth(1).unwrap().parse::<String>().unwrap();
    
    if mode == "compute" {
        let theta = FixedPoint::new(std::env::args().nth(2).unwrap().parse::<f64>().unwrap());
        let iters = std::env::args().nth(3).unwrap().parse::<usize>().unwrap();

        let ret = cordic(theta, iters);
        println!("cos {} == {}\nsin {} == {}", theta, ret[0], theta, ret[1]);
    } else if mode == "bench" {
        // NOTE: Output is a CSV file that I will open in Excel
        println!("Theta, CORDIC Cosine, Standard Cosine, Cosine Error, CORDIC Sine, Standard Sine, Sine Error");
        for i in 0..314 {
            let theta = i as f64 / 100.0;

            let cordic_val = cordic(FixedPoint::new(theta), 100);
            let cos_val = (theta).cos();
            let sin_val = (theta).sin();

            // NOTE: Pulling from `val` directly is a hack and shouldn't be allowed
            println!(
                "{},{},{},{},{},{}", 
                cordic_val[0], cos_val, (cordic_val[0].val - cos_val).abs(),
                cordic_val[1], sin_val, (cordic_val[1].val - sin_val).abs()
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close_enough(a: FixedPoint, b: FixedPoint) -> bool {
        (a - b) < FixedPoint::new(0.01)
    }

    #[test]
    fn basic() {
        for i in 0..157 {
            // pi/2
            let ret = cordic(FixedPoint::new(i as f64 / 100.0), 1000);
            let cos = FixedPoint::new((i as f64 / 100.0).cos());
            let sin = FixedPoint::new((i as f64 / 100.0).sin());

            println!(
                "Theta == {}\t{} vs {}\t{} vs {}",
                (i as f64 / 100.0),
                ret[0],
                cos,
                ret[1],
                sin
            );

            assert![close_enough(ret[0], cos)];
            assert![close_enough(ret[1], sin)];
        }
    }
}
