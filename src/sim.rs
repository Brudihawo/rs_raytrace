use serde::{Deserialize, Serialize};

use crate::ray::{RayTemplate, Ray};

use crate::util::*;
use crate::boundary::BoundaryType;

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
pub struct SimTemplate {
    mode: Mode,
    ray: RayTemplate,
}

pub struct RaySim {
    mode: Mode,
    ray: Ray,
}

impl RaySim {
    pub fn from_template(template: SimTemplate) -> RaySim {
        RaySim {
            mode: template.mode,
            ray: Ray::from_template(template.ray),
        }
    }

    pub fn sim(&mut self) {
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
                    max_angle = deg_to_rad(max_angle);
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
                height: radius, opt_idx, ..
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
