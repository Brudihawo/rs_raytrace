use std::f64::consts::FRAC_PI_2;

use serde::{Deserialize, Serialize};
use crate::boundary::BoundaryType;

use crate::util::*;

#[derive(Serialize, Deserialize)]
pub struct RayTemplate {
    x: f64,
    y: f64,
    angle: f64,
    o_idx: f64,
    boundaries: Vec<BoundaryType>,
}

#[derive(Serialize, Deserialize)]
pub struct Ray {
    // Original State
    pub o_x: f64,
    pub o_y: f64,
    pub o_angle: f64,
    pub o_n: f64,
    // Current State
    pub x: f64,
    pub y: f64,
    pub angle: f64,
    pub cur_idx: f64,
    pub boundary: usize,
    pub boundaries: Vec<BoundaryType>,
}

impl Ray {
    pub fn from_template(template: RayTemplate) -> Ray {
        Ray {
            o_x: template.x,
            o_y: template.y,
            o_n: template.o_idx,
            o_angle: template.angle,

            x: template.x,
            y: template.y,
            cur_idx: template.o_idx,
            angle: 0.0,

            boundary: 0,
            boundaries: template.boundaries
        }
    }
    fn check_total_reflection(angle: f64, cur_n: f64, next_n: f64) -> bool {
        if next_n < cur_n {
            let alpha_g = (next_n / cur_n).asin().abs();
            return angle.abs() > alpha_g.abs();
        }
        false
    }
    pub fn reset(&mut self, angle: f64, y: f64) {
        self.x = self.o_x;
        self.y = y;
        self.cur_idx = self.o_n;
        self.angle = angle;
        self.o_angle = angle;
        self.boundary = 0;
    }

    pub fn print_state(&self, verbose: bool, id: usize) {
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

    pub fn next_boundary(&mut self) -> bool {
        if self.boundary == self.boundaries.len() {
            return false;
        }
        let bound = self.boundary;
        self.boundary += 1;
        match &self.boundaries[bound] {
            BoundaryType::Line {
                opt_idx,
                midpoint,
                height,
            } => {
                let delta_x = midpoint - self.x;
                let delta_y = delta_x * self.angle.tan();

                // Check if out of bounds
                if (self.y + delta_y).abs() > *height {
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

