use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub enum BoundaryType {
    Line {
        opt_idx: f64,
        midpoint: f64,
        height: f64, // maximum distance from optical axis
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
