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

// Fixed-Point Arithmetic
//
// The only reason this is implemented is to showcase the
// efficiency of bitshifts-as-multiplication. All functions
// done through here can be done through f64 normally
pub struct FixedPoint {
    val: f64
}

impl FixedPoint {
    fn new(val: i64) -> Self {
	Self {
	    val: val as f64
	}
    }
}

impl Add for FixedPoint {
    type Output = Self;
    fn add(self, other: Self) -> Self {
	Self { val: self.val + other.val }
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

fn main() {
    // Pull parameters from string, should be called as either
    // ./cordic-rs [theta] [iters]
    // or
    // cargo run [theta] [iters]
    let theta = std::env::args().nth(1).unwrap().parse::<f32>().unwrap();
    let iters = std::env::args().nth(2).unwrap().parse::<u64>().unwrap();

    // CORDIC (for trig functions, at least) does require some
    // compile time constants. However, this is far more space
    // efficient than naively storing sine itself.
    
    // atan(2^-x)
    let angles = (0_i32..28_i32).map(|x| {
	(2_f32.powi(-1 * x) as f32).atan()
    }).collect::<Vec<f32>>();

    // cumprod(1 / sqrt(1 + 2^-2y))
    // NOTE The cumulative product is done by re-calculating and multiplying
    // all elements of the vector together with fold()
    let kvalues = (0..23).map(|x| {
	(0_i32..23_i32).take(x + 1).map(|y| {
	    1.0_f32 / (1.0_f32 + 2_f32.powi(-2 * y)).sqrt().abs()
	}).fold(1.0, |x, y| x*y)
    }).collect::<Vec<f32>>();

    let mut poweroftwo = 1.0;
    let mut angle = angles[0];
    let mut v = [1.0, 0.0]; // Initialize as cos = 1, sine = 0
    let mut cur_theta = theta;
    for i in 1..iters-1 {
	let sigma = if cur_theta < 0.0 {
	    -1.0
	} else {
	    1.0
	};

	let factor = sigma * poweroftwo;
	let matrix = [
	    [1.0, -factor],
	    [factor, 1.0]
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
	poweroftwo = poweroftwo / 2.0;

	if i + 2 > angles.len() as u64 {
	    angle = angle / 2.0;
	} else {
	    angle = angles[i as usize + 2];
	}
    }

    // Scale vector back such that magnitude is 1
    // NOTE: This can be done either by keeping track of the
    // initial values or performing a square root
    v = [
	v[0] * kvalues[iters as usize - 1],
	v[1] * kvalues[iters as usize - 1]
    ];
    v = [
	v[0] / (v[0]*v[0] + v[1]*v[1]).sqrt(),
	v[1] / (v[0]*v[0] + v[1]*v[1]).sqrt(),
    ];
    println!("cos {} == {}", theta, v[0]);
    println!("sin {} == {}", theta, v[1]);
}
