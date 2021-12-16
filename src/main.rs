use std::f64::consts::{FRAC_1_PI, PI};

fn deg_to_rad(deg: f64) -> f64 {
    PI / 180.0f64 * deg
}

fn rad_to_deg(rad: f64) -> f64 {
    180.0f64 * FRAC_1_PI * rad
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
        angle_radius: f64,
    },
}

struct Ray {
    x: f64,
    y: f64,
    angle: f64,
    cur_n: f64,
    boundaries: Vec<BoundaryType>,
}

impl Ray {
    fn next_boundary(&mut self) {
        for boundary in self.boundaries {
            match boundary {
                BoundaryType::Line {
                    opt_idx,
                    midpoint,
                    radius,
                } => {
                    let delta_x = midpoint - self.x;
                    let delta_y = delta_x * self.angle.tan();
                    self.x += delta_x;
                    self.y += delta_y;
                }
                BoundaryType::Spherical {
                    opt_idx,
                    midpoint,
                    radius,
                    angle_radius,
                } => {}
            }
        }
    }
}

fn main() {
    let boundaries = Vec::from([
        BoundaryType::Spherical {
            opt_idx: 2.1,
            midpoint: 10.0,
            radius: 5.0,
            angle_radius: deg_to_rad(50.0),
        },
        BoundaryType::Line {
            opt_idx: 1.0,
            midpoint: 15.0,
            radius: 10.0,
        },
    ]);

    let ray: Ray = Ray {
        x: 0.0,
        y: 0.0,
        angle: deg_to_rad(10.0),
        cur_n: 1.0,
        boundaries,
    };
}
