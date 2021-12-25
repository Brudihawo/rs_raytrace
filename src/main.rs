use serde::{Deserialize, Serialize};
use serde_json;
use std::f64::consts::{FRAC_1_PI, FRAC_PI_2, PI};
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

fn conic(r: f64, kappa: f64, radius: f64) -> f64 {
    r.powi(2) / (radius * (1.0 - (1.0 + kappa) * r.powi(2) / radius.powi(2)))
}

#[derive(Serialize, Deserialize)]
enum BoundaryType {
    Line {
        opt_idx: f64,
        midpoint: f64,
        radius: f64, // maximum distance from optical axis
    },
    Spherical {
        opt_idx: f64,
        midpoint: f64,
        radius: f64,
        height: f64, // maximum distance from optical axis
    },
    Conic {
        opt_idx: f64,
        midpoint: f64,
        radius: f64,
        conic_param: f64,
        height: f64, // maximum distance from optical axis
    },
}

#[derive(Serialize, Deserialize)]

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
    cur_idx: f64,
    boundary: usize,
    boundaries: Vec<BoundaryType>,
}

impl Ray {
    fn check_total_reflection(angle: f64, cur_n: f64, next_n: f64) -> bool {
        if next_n < cur_n {
            let alpha_g = (next_n / cur_n).asin().abs();
            return angle.abs() > alpha_g.abs();
        }
        false
    }
    fn new(x: f64, y: f64, angle: f64, n: f64, boundaries: Vec<BoundaryType>) -> Ray {
        Ray {
            x,
            y,
            angle,
            cur_idx: n,
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
        self.cur_idx = self.o_n;
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
                self.cur_idx
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
                    if Ray::check_total_reflection(self.angle, self.cur_idx, *opt_idx) {
                        return false;
                    }

                    // compute geometric refraction
                    self.angle = (self.cur_idx / opt_idx * self.angle.sin()).asin();
                    self.cur_idx = *opt_idx;

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
                    if Ray::check_total_reflection(alpha, self.cur_idx, *opt_idx) {
                        return false;
                    }

                    // compute geometric refraction
                    self.angle = inter_angle - (self.cur_idx / opt_idx * alpha.sin()).asin();
                    self.cur_idx = *opt_idx;

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
                    if Ray::check_total_reflection(alpha, self.cur_idx, *opt_idx) {
                        return false;
                    }

                    // compute geometric refraction
                    self.angle = inter_angle - (self.cur_idx / opt_idx * alpha.sin()).asin();
                    self.cur_idx = *opt_idx;

                    if self.angle.is_nan() {
                        return false;
                    }

                    return true;
                }
            }
            BoundaryType::Conic {
                opt_idx,
                midpoint,
                radius,
                conic_param,
                height,
            } => {
                let x0 = self.x;
                let y0 = self.y;
                let phi = self.angle;

                let c0 = -(1.0 + conic_param + phi.tan().powi(2));
                let c1 =
                    2.0 * (radius - (1.0 + conic_param) * (x0 - midpoint) - y0 * self.angle.tan());
                let c2 = 2.0 * radius * (x0 - midpoint)
                    - (1.0 + conic_param) * (x0 - midpoint).powi(2)
                    - y0.powi(2);

                let p = c1 / c0;
                let q = c2 / c0;
                let det = (p / 2.0).powi(2) - q;

                // No solution / no intersection
                if det < 0.0 {
                    return false;
                }

                let mut a: f64;
                if *radius > 0.0 {
                    if *conic_param >= -1.0 {
                        a = -p / 2.0 - det.sqrt();
                    } else {
                        a = -p / 2.0 + det.sqrt();
                    }
                } else {
                    a = -p / 2.0 + det.sqrt();
                }

                // Handle intersection on optical axis and ray angle 0
                if c0 == 0.0 {
                    a = -c2 / c1;
                }

                self.x += a;
                self.y += a * self.angle.tan();

                // Out of bounds check
                if self.y.abs() > *height {
                    return false;
                }

                let inter_angle: f64; // angle of boundary at intersection
                let alpha: f64; // angle of incidence

                let diff = (radius - (1.0 + conic_param) * (self.x - midpoint))
                    / (2.0 * (self.x - midpoint) * radius
                        - (self.x - midpoint).powi(2) * (1.0 + conic_param))
                        .sqrt();

                if self.y < 0.0 {
                    if *radius > 0.0 {
                        inter_angle = -diff.atan();
                        alpha = -FRAC_PI_2 - inter_angle + self.angle;
                        self.angle =
                            FRAC_PI_2 + inter_angle + (self.cur_idx / opt_idx * alpha.sin()).asin();
                    } else {
                        inter_angle = diff.atan();
                        alpha = -FRAC_PI_2 - inter_angle - self.angle;
                        self.angle = -FRAC_PI_2
                            - inter_angle
                            - (self.cur_idx / opt_idx * alpha.sin()).asin();
                    }
                    // incident angle
                } else if self.y > 0.0 {
                    if *radius > 0.0 {
                        inter_angle = diff.atan();
                        alpha = FRAC_PI_2 - inter_angle + self.angle;
                        self.angle = -FRAC_PI_2
                            + inter_angle
                            + (self.cur_idx / opt_idx * alpha.sin()).asin();
                    } else {
                        inter_angle = (-diff).atan();
                        alpha = FRAC_PI_2 - inter_angle - self.angle;
                        self.angle =
                            FRAC_PI_2 - inter_angle - (self.cur_idx / opt_idx * alpha.sin()).asin();
                    }
                } else {
                    // handle intersection at y=0
                    alpha = self.angle; // for debug purposes
                    self.angle = (self.cur_idx / opt_idx * self.angle.sin()).asin();
                    inter_angle = FRAC_PI_2; // for debug purposes
                }
                // check for total reflection (indirectly)
                if self.angle.is_nan() {
                    return false;
                }
                self.cur_idx = *opt_idx;

                true
            }
        }
    }
}

#[derive(Serialize, Deserialize)]
enum Mode {
    RayFan {
        n_rays: usize,
        min_angle: f64,
        max_angle: f64,
        in_deg: bool,
    },
    RayArray {
        n_rays: usize,
        min_height: f64,
        max_height: f64,
    },
}

#[derive(Serialize, Deserialize)]
struct RaySim {
    mode: Mode,
    ray: Ray,
}

impl RaySim {
    fn sim(&mut self) {
        match self.mode {
            Mode::RayArray {
                n_rays,
                min_height,
                max_height,
            } => {
                for idx in 0..n_rays + 1 {
                    let cur_height = lerp(min_height, max_height, idx as f64 / n_rays as f64);
                    self.ray.reset(self.ray.o_angle, cur_height);
                    self.ray.print_state(false, idx);
                    while self.ray.next_boundary() {
                        self.ray.print_state(false, idx);
                    }
                    println!("");
                }
            }
            Mode::RayFan {
                n_rays,
                mut min_angle,
                mut max_angle,
                in_deg,
            } => {
                if in_deg {
                    min_angle = deg_to_rad(min_angle);
                    max_angle = deg_to_rad(min_angle);
                }
                for idx in 0..n_rays + 1 {
                    // convert angle
                    let cur_angle = lerp(min_angle, max_angle, idx as f64 / n_rays as f64);
                    self.ray.reset(cur_angle, self.ray.o_y);
                    self.ray.print_state(false, idx);
                    while self.ray.next_boundary() {
                        self.ray.print_state(false, idx);
                    }
                    println!("");
                }
            }
        }
    }
}

fn read_infile(fname: String) -> Result<(RaySim, Mode), ()> {
    Err(())
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
            BoundaryType::Conic {
                radius,
                opt_idx,
                conic_param,
                height,
                ..
            } => {
                if radius == 0.0 || opt_idx < 0.0 {
                    eprintln!(
                        "Invalid Parameters for Conic Boundary at index {}: Invalid values",
                        idx
                    );
                    return false;
                }
                // maximum reachable radius is larger than lens parametrisation allows
                if height.powi(2) > radius.powi(2) / (1.0 + conic_param) && conic_param > -1.0 {
                    eprintln!(
                        "Invalid Parameters for Conic Boundary at index {}.\n\
                               Height exceeds maximum value allowed by lens parametrisation.\n\
                               H < sqrt(R^2 / (1 + conic_param)) (Here: H < {})",
                        idx,
                        (radius.powi(2) / (1.0 + conic_param)).sqrt()
                    );
                    return false;
                }
            }
        }
    }
    true
}

fn ray_fan(mut ray: Ray, n_rays: usize, min_angle: f64, max_angle: f64, in_deg: bool) {
    for idx in 0..(n_rays + 1) {
        ray.reset(
            lerp(
                if in_deg {
                    deg_to_rad(min_angle)
                } else {
                    min_angle
                },
                if in_deg {
                    deg_to_rad(max_angle)
                } else {
                    max_angle
                },
                idx as f64 / n_rays as f64,
            ),
            0.0,
        );
        ray.print_state(false, idx);
        while ray.next_boundary() {
            ray.print_state(false, idx);
        }
        println!("");
    }
}

fn ray_array(mut ray: Ray, n_rays: usize, min_height: f64, max_height: f64) {
    for idx in 0..(n_rays + 1) {
        ray.reset(
            0.0,
            lerp(min_height, max_height, idx as f64 / n_rays as f64),
        );
        ray.print_state(false, idx);
        while ray.next_boundary() {
            ray.print_state(false, idx);
        }
        println!("");
    }
}

fn main() {
    let boundaries1 = Vec::from([
        BoundaryType::Conic {
            opt_idx: 1.5,
            midpoint: 15.0,
            radius: 10.0,
            conic_param: -0.9,
            height: 10.0,
        },
        BoundaryType::Conic {
            opt_idx: 1.0,
            midpoint: 20.0,
            radius: -10.0,
            conic_param: -0.9,
            height: 10.0,
        },
        BoundaryType::Line {
            opt_idx: 2.0,
            midpoint: 43.0,
            radius: 10.0,
        },
    ]);

    let boundaries2 = Vec::from([
        BoundaryType::Conic {
            opt_idx: 1.5,
            midpoint: 10.0,
            radius: 10.0,
            conic_param: -0.9,
            height: 10.0,
        },
        BoundaryType::Line {
            opt_idx: 1.0,
            midpoint: 20.0,
            radius: 10.0,
        },
        BoundaryType::Line {
            opt_idx: 2.0,
            midpoint: 28.0,
            radius: 20.0,
        },
    ]);

    let boundaries3 = Vec::from([
        BoundaryType::Line {
            opt_idx: 1.5,
            midpoint: 15.0,
            radius: 10.0,
        },
        BoundaryType::Conic {
            opt_idx: 1.0,
            midpoint: 20.0,
            radius: -10.0,
            conic_param: -0.9,
            height: 10.0,
        },
        BoundaryType::Line {
            opt_idx: 2.0,
            midpoint: 28.0,
            radius: 20.0,
        },
    ]);

    let mut sim = RaySim {
        mode: Mode::RayArray {
            n_rays: 20,
            max_height: 4.0,
            min_height: -4.0,
        },
        ray: Ray::new(0.0, 0.0, 0.0, 1.0, boundaries1),
    };

    sim.sim();

    {
        let mut file = File::create("boundaries.json").unwrap();
        let mut sim_json_writer = serde_json::to_writer_pretty(&file, &sim).unwrap();
        file.write(b"\n");
    }
}

