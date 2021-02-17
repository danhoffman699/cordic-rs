// CORDIC Demonstration in Rust
//
// If you have any questions, feel free to email me at djh4@illinois.edu

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
	(0_i32..23_i32).take(x).map(|y| {
	    1.0_f32 / (1.0_f32 + 2_f32.powi(-2 * y)).sqrt().abs()
	}).fold(1.0, |x, y| x*y)
    }).collect::<Vec<f32>>();
    
    println!("angles: {:#?}", &angles);
    println!("kvalues: {:#?}", &kvalues);

    println!("theta: {}, iters: {}", theta, iters);

    let mut poweroftwo = 1.0;
    let mut angle = angles[1];
    let mut v = [1.0, 0.0]; // Initialize as cos = 1, sine = 0
    let mut cur_theta = theta;
    for i in 0..iters-1 {
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
	
	v = [
	    matrix[0][0] * v[0] + matrix[0][1] * v[1],
	    matrix[1][0] * v[0] + matrix[1][1] * v[1]
	];

	cur_theta = theta - sigma * angle;
	poweroftwo = poweroftwo / 2.0;

	if i + 2 > angles.len() as u64 {
	    angle = angle / 2.0;
	} else {
	    angle = angles[i as usize + 2];
	}
    }
}
