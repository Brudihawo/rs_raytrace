use std::f64::consts::{FRAC_1_PI, PI};
use std::fs::File;
use std::io::prelude::*;

fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0f64
}

fn rad_to_deg(rad: f64) -> f64 {
    180.0f64 * FRAC_1_PI * rad
}

fn lerp(x0: f64, x1: f64, t: f64) -> f64 {
    return x0 * (1.0 - t) + x1 * t;
}

enum BoundaryType {
    Line {
        opt_idx: f64,
        midpoint: f64,
        radius: f64,
    },
    Spherical {
        opt_idx: f64,
        midpoint: f64,
        radius: f64,
        height: f64,
    },
}

struct Ray {
    // Original State
    o_x: f64,
    o_y: f64,
    o_angle: f64,
    o_n: f64,
    // Current State
    x: f64,
    y: f64,
    angle: f64,
    cur_n: f64,
    boundary: usize,
    boundaries: Vec<BoundaryType>,
}

impl Ray {
    fn check_total_reflection(angle: f64, cur_n: f64, next_n: f64) -> bool {
        if next_n < cur_n {
            let alpha_g = (cur_n / next_n).asin().abs();
            return angle.abs() > alpha_g.abs();
        }
        false
    }
    fn new(x: f64, y: f64, angle: f64, n: f64, boundaries: Vec<BoundaryType>) -> Ray {
        Ray {
            x,
            y,
            angle,
            cur_n: n,
            boundary: 0,
            boundaries,

            o_x: x,
            o_y: y,
            o_angle: angle,
            o_n: n,
        }
    }

    fn reset(&mut self, angle: f64, y: f64) {
        self.x = self.o_x;
        self.y = y;
        self.cur_n = self.o_n;
        self.angle = angle;
        self.o_angle = angle;
        self.boundary = 0;
    }

    fn print_state(&self, verbose: bool, id: usize) {
        if verbose {
            println!(
                "Hit Boundary {} at ({}, {}) | new angle: {} deg, new n = {}",
                self.boundary,
                self.x,
                self.y,
                rad_to_deg(self.angle),
                self.cur_n
            );
        } else {
            println!("{}, {}, {}, {}", id, self.x, self.y, rad_to_deg(self.angle));
        }
    }

    fn next_boundary(&mut self) -> bool {
        if self.boundary == self.boundaries.len() {
            return false;
        }
        let bound = self.boundary;
        self.boundary += 1;
        match &self.boundaries[bound] {
            BoundaryType::Line {
                opt_idx,
                midpoint,
                radius,
            } => {
                let delta_x = midpoint - self.x;
                let delta_y = delta_x * self.angle.tan();

                // Check if out of bounds
                if (self.y + delta_y).abs() > *radius {
                    return false;
                } else {
                    self.x += delta_x;
                    self.y += delta_y;

                    // test total reflection
                    if Ray::check_total_reflection(self.angle, self.cur_n, *opt_idx) {
                        // println!("Total Reflection");
                        return false;
                    }

                    // compute geometric refraction
                    self.angle = (self.cur_n / opt_idx * self.angle.sin()).asin();
                    self.cur_n = *opt_idx;

                    return true;
                }
            }

            BoundaryType::Spherical {
                opt_idx,
                midpoint,
                radius,
                height,
            } => {
                let tmp_val = self.angle.cos().abs() / radius
                    * (self.angle.tan() * (midpoint - self.x) + self.y);
                if tmp_val.abs() > 1.0 {
                    return false;
                }

                let inter_angle = -tmp_val.asin() + self.angle;

                if *radius > 0.0 {
                    self.x = midpoint - radius * inter_angle.cos();
                    self.y = -radius * inter_angle.sin();

                    // out of bounds check
                    if self.y.abs() > *height {
                        return false;
                    }

                    // incident angle
                    let alpha = inter_angle - self.angle;
                    if Ray::check_total_reflection(alpha, self.cur_n, *opt_idx) {
                        return false;
                    }

                    // compute geometric refraction
                    self.angle = inter_angle - (self.cur_n / opt_idx * alpha.sin()).asin();
                    self.cur_n = *opt_idx;

                    if self.angle.is_nan() {
                        return false;
                    }

                    return true;
                } else {
                    self.x = midpoint - radius * inter_angle.cos();
                    self.y = -radius * inter_angle.sin();

                    // out of bounds check
                    if self.y.abs() > *height {
                        return false;
                    }

                    // incident angle
                    let alpha = inter_angle - self.angle;
                    if Ray::check_total_reflection(alpha, self.cur_n, *opt_idx) {
                        return false;
                    }

                    // compute geometric refraction
                    self.angle = inter_angle - (self.cur_n / opt_idx * alpha.sin()).asin();
                    self.cur_n = *opt_idx;

                    if self.angle.is_nan() {
                        return false;
                    }

                    return true;
                }
            }
        }
    }
}

fn dump_boundaries(boundaries: &Vec<BoundaryType>) -> Result<(), std::io::Error> {
    let fname = "boundaries.json";
    let mut file = File::create(fname)?;
    file.write(b"{\n")?;
    for (i, boundary) in boundaries.iter().enumerate() {
        file.write(format!("  \"{}\": {{\n", i).as_bytes())?;
        match boundary {
            BoundaryType::Line {
                opt_idx,
                midpoint,
                radius,
            } => {
                file.write(
                    format!(
                        "    \"type\": \"line\",\n    \
                             \"opt_idx\": {},\n    \
                             \"midpoint\": {},\n    \
                             \"radius\": {}\n",
                        opt_idx, midpoint, radius
                    )
                    .as_bytes(),
                )?;
            }
            BoundaryType::Spherical {
                opt_idx,
                midpoint,
                radius,
                height,
            } => {
                file.write(
                    format!(
                        "    \"type\": \"spherical\",\n    \
                             \"opt_idx\": {},\n    \
                             \"midpoint\": {},\n    \
                             \"radius\": {},\n    \
                             \"height\": {}\n",
                        opt_idx, midpoint, radius, height
                    )
                    .as_bytes(),
                )?;
            }
        }
        if i < boundaries.len() - 1 {
            file.write(b"  },\n")?;
        } else {
            file.write(b"  }\n")?;
        }
    }
    file.write(b"}\n")?;
    Ok(())
}

fn validate_boundaries(boundaries: &Vec<BoundaryType>) -> bool {
    for (idx, boundary) in boundaries.iter().enumerate() {
        match *boundary {
            BoundaryType::Spherical {
                radius,
                opt_idx,
                height,
                ..
            } => {
                if radius == 0.0 || height < 0.0 || opt_idx < 0.0 {
                    eprintln!("Invalid Parameters for Spherical Boundary at index {}", idx);
                    return false;
                }
            }
            BoundaryType::Line {
                radius, opt_idx, ..
            } => {
                if radius == 0.0 || opt_idx < 0.0 {
                    eprintln!("Invalid Parameters for Line Boundary at index {}", idx);
                    return false;
                }
            }
        }
    }
    for (idx, pair) in boundaries.windows(2).enumerate() {
        let _b0_rightmost: f64;
        let _b1_leftmost: f64;
        let b0_idx: f64;
        let b1_idx: f64;

        match pair[0] {
            BoundaryType::Line {
                opt_idx, midpoint, ..
            } => {
                _b0_rightmost = midpoint;
                b0_idx = opt_idx;
            }
            BoundaryType::Spherical {
                opt_idx,
                midpoint,
                radius,
                ..
            } => {
                _b0_rightmost = midpoint - radius;
                b0_idx = opt_idx;
            }
        }

        match pair[1] {
            BoundaryType::Line {
                opt_idx, midpoint, ..
            } => {
                _b1_leftmost = midpoint;
                b1_idx = opt_idx;
            }
            BoundaryType::Spherical {
                opt_idx,
                midpoint,
                radius,
                height,
            } => {
                if radius < 0.0 {
                    _b1_leftmost = midpoint + (radius.powi(2) - height.powi(2)).sqrt()
                } else {
                    _b1_leftmost = midpoint - radius;
                }
                b1_idx = opt_idx;
            }
        }

        if b0_idx == b1_idx {
            eprintln!(
                "Optical indices are the same at boundaries {} and {}!",
                idx,
                idx + 1
            );
            return false;
        }
    }
    false
}

fn main() {
    let boundaries = Vec::from([
        BoundaryType::Line {
            opt_idx: 2.0,
            midpoint: 5.0,
            radius: 10.0,
        },
        BoundaryType::Line {
            opt_idx: 1.0,
            midpoint: 10.0,
            radius: 10.0,
        },
        BoundaryType::Spherical {
            opt_idx: 2.0,
            midpoint: 30.0,
            radius: 15.0,
            height: 8.0,
        },
        // BoundaryType::Line {
        //     opt_idx: 2.0,
        //     midpoint: 20.0,
        //     radius: 10.0,
        // },
        BoundaryType::Spherical {
            opt_idx: 1.1,
            midpoint: 5.0,
            radius: -15.0,
            height: 8.0,
        },
        // BoundaryType::Spherical {
        //     opt_idx: 1.0,
        //     midpoint: 32.0,
        //     radius: 15.0,
        //     height: 10.0,
        // },
        BoundaryType::Line {
            opt_idx: 2.0,
            midpoint: 32.0,
            radius: 10.0,
        },
    ]);

    dump_boundaries(&boundaries).unwrap();

    let mut ray: Ray = Ray::new(0.0, 0.0, deg_to_rad(20.0), 1.0, boundaries);

    const N_RAYS: usize = 20;
    for idx in 0..(N_RAYS + 1) {
        // ray.reset(
        //     lerp(
        //         deg_to_rad(-30.0),
        //         deg_to_rad(30.0),
        //         idx as f64 / N_RAYS as f64,
        //     ),
        //     0.0,
        // );
        ray.reset(0.0, lerp(-5.0, 5.0, idx as f64 / N_RAYS as f64));
        ray.print_state(false, idx);
        while ray.next_boundary() {
            ray.print_state(false, idx);
        }
        println!("");
    }
}
