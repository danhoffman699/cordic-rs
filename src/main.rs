// CORDIC Demonstration in Rust
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

use std::ops::{ Add, Sub, Mul, Div };
use std::cmp::{ Ordering, PartialEq, PartialOrd };
use std::fmt;

// Fixed-Point Arithmetic
//
// The only reason this is implemented is to showcase the
// efficiency of bitshifts-as-multiplication. All functions
// done through here can be done through f64 normally
pub struct FixedPoint {
    val: fixed::types::I4F124
}

impl FixedPoint {
    fn new(val: f64) -> Self {
	Self {
	    val: fixed::types::I4F124::from_num(val)
	}
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
	Self {val: self.val + other.val }
    }
}

impl Sub for FixedPoint {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
	Self { val: self.val - other.val }
    }
}

impl Mul for FixedPoint {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
	Self { val: self.val * other.val }
    }
}

impl Div for FixedPoint {
    type Output = Self;
    fn div(self, other: Self) -> Self {
	Self { val: self.val / other.val }
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

fn cordic(theta: FixedPoint, iters: u64) -> [FixedPoint; 2] {
    // CORDIC (for trig functions, at least) does require some
    // compile time constants. However, this is far more space
    // efficient than naively storing sine itself.
    
    // atan(2^-x)
    let angles = (0_i32..28_i32).map(|x| {
	FixedPoint::new(
	    (2_f64.powi(-1 * x) as f64).atan()
	)
    }).collect::<Vec<FixedPoint>>();

    // cumprod(1 / sqrt(1 + 2^-2y))
    // NOTE The cumulative product is done by re-calculating and multiplying
    // all elements of the vector together with fold()
    let kvalues = (0..23).map(|x| {
	FixedPoint::new(
	    (0_i32..23_i32).take(x + 1).map(|y| {
		1.0_f64 / (1.0_f64 + 2_f64.powi(-2 * y)).sqrt().abs()
	    }).fold(1.0, |x, y| x*y)
	)
    }).collect::<Vec<FixedPoint>>();

    let mut poweroftwo = FixedPoint::new(1.0);
    let mut angle = angles[0];
    let mut v = [FixedPoint::new(1.0), FixedPoint::new(0.0)]; // Initialize as cos = 1, sine = 0
    let mut cur_theta = theta;
    for i in 1..iters-1 {
	let sigma = FixedPoint::new(
	    if cur_theta < FixedPoint::new(0.0) {
		-1.0
	    } else {
		1.0
	    }
	);

	let factor = sigma * poweroftwo;
	let matrix = [
	    [FixedPoint::new(1.0), FixedPoint::new(-1.0) * factor],
	    [factor, FixedPoint::new(1.0)]
	];
	
	// v = R * v
	// NOTE: Matrix is always of the form
	// [ 1.0, -factor; factor, 1.0 ]
	//
	// The statements x *= 2 and x <<= 1 are equivalent
	v = [
	    matrix[0][0] * v[0] + matrix[0][1] * v[1],
	    matrix[1][0] * v[0] + matrix[1][1] * v[1]
	];

	cur_theta = cur_theta - sigma * angle;
	poweroftwo = poweroftwo / FixedPoint::new(2.0);

	if i + 2 > angles.len() as u64 {
	    angle = angle / FixedPoint::new(2.0);
	} else {
	    angle = angles[i as usize + 2];
	}
    }

    // Scale vector back such that magnitude is 1
    // NOTE: This can be done either by keeping track of the
    // initial values or performing a square root. If the machine
    // is slow enough that CORDIC is practical (i.e. expensive
    // hardware multiplication), then it is too slow for square roots
    // and divisions
    v = [
	v[0] * kvalues[iters as usize - 1],
	v[1] * kvalues[iters as usize - 1]
    ];

    v
}

fn main() {
    // Pull parameters from string, should be called as either
    // ./cordic-rs [theta] [iters]
    // or
    // cargo run [theta] [iters]
    let theta = FixedPoint::new(std::env::args().nth(1).unwrap().parse::<f64>().unwrap());
    let iters = std::env::args().nth(2).unwrap().parse::<u64>().unwrap();

    cordic(theta, iters);
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn basic() {
	let ret = cordic(FixedPoint::new(0.0), 10);

	assert![ret[0] > FixedPoint::new(0.8) && ret[0] < FixedPoint::new(1.2)];
	assert![ret[1] > FixedPoint::new(-0.2) && ret[1] < FixedPoint::new(0.2)];
    }
}
