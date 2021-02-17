// CORDIC Demonstration in Rust
//
// If you have any questions, feel free to email me at djh4@illinois.edu

fn main() {
    // CORDIC (for trig functions, at least) does require some
    // compile time constants. However, this is far more space
    // efficient than naively storing sine itself.
    
    let angles = (0_i32..28_i32).map(|x| {
	(2_f32.powi(-1 * x) as f32).atan()
    }).collect::<Vec<f32>>();

    let kvalues = (0..23).map(|x| {
	(0_i32..23_i32).take(x).map(|y| {
	    1.0_f32 / (1.0_f32 + 2_f32.powi(-2 * y)).sqrt().abs()
	}).fold(1.0, |x, y| x*y)
    }).collect::<Vec<f32>>();
    
    println!("angles: {:#?}", &angles);
    println!("kvalues: {:#?}", &kvalues);

    
}
