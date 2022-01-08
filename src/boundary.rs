use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub enum BoundaryType {
    Line {
        opt_idx: f64,      // optical index after boundary
        x: f64,            // coordinate in x
        y: f64,            // coordinate in y
        height: f64,       // maximum distance from point
        angle: f64,        // angle to vertical line
    },
    Spherical {
        opt_idx: f64,      // optical index after boundary
        x: f64,            // coordinate in x
        // y: f64,            // coordinate in y
        radius: f64,       // circular radius
        height: f64,       // maximum distance from optical axis
    },
    Conic {
        opt_idx: f64,
        midpoint: f64,
        radius: f64,
        conic_param: f64,
        height: f64,       // maximum distance from optical axis
    },
}
