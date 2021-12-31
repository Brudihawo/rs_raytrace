use std::f64::consts::{FRAC_1_PI, PI};

pub fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0f64
}

pub fn rad_to_deg(rad: f64) -> f64 {
    180.0f64 * FRAC_1_PI * rad
}

pub fn lerp(x0: f64, x1: f64, t: f64) -> f64 {
    return x0 * (1.0 - t) + x1 * t;
}

pub fn exit(message: String, errcode: i32) {
    eprintln!("{}", message);
    std::process::exit(errcode);
}
