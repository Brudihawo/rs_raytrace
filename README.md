# Ray Tracing Demo for Geometric Optics in Rust

## Roadmap
1. Implement mirrors
2. Ray Energy Content via reflection / transmission
3. Optimisaiton of Key features

## Quickstart
```commandline
  $ git clone https://github.com/brudihawo/rs_raytrace
  $ cd rs_raytrace
  $ cargo build
```

To execute an example simulation and plot (in project root):
```commandline
  $ ./scripts/exec_plot.sh example_infile.json
```

## Features
- Mathematically accurate raytracing
- Lens Shapes:
  - Vertical Line
  - Spherical
  - Conic
- Plotting Visualisation using Python + Matplotlib
- Read simulation configuration from file

# Infile Syntax
Simulation config is provided as JSON. The overall structure is:
```json
{
  <mode>,
  <ray>
}
```

`mode` can be either
```json
"RayFan": {
  "n_rays": <usize>,
  "min_angle": <float>,
  "max_angle": <float>,
  "in_deg": bool
}
```
or
```json
"RayArray": {
  "n_rays": <usize>,
  "min_height": <float>,
  "max_height": <float>,
}
```
and specifies the way the Ray parameters are varied during the simulation process.
`RayFan` casts a fan of rays with varying angles from the same position. `RayArray` casts
an array of rays with the same angle from the same x coordinate, but with varying y
coordinate (distance to the optical axis).

`<ray>` is a specification of the Rays general properties:
```json
"ray": {
  "x": <float>,                  // initial x coordinate
  "y": <float>,                  // initial y coordinate
  "angle": <float>,              // initial ray angle
  "o_idx": <float>,              // initial refractive index
  "boundaries": [<BoundaryType>] // boundaries
}
```
`boundaries` is an array of boundaries along the optical axis, specifying the shape of
the boundaries in sequential order and changes in refractive index at those borders.
Currently, the following boundary types are implemented:
```json
"Line": {
  "midpoint": <float>,
  "height": <float>,
  "opt_idx": <float>,
},

"Spherical": {
  "midpoint": <float>,
  "radius": <float>,
  "height": <float>,
  "opt_idx": <float>,
},

"Conic": {
  "midpoint": <float>,
  "radius": <float>,
  "conic_param": <float>,
  "height": <float>,
  "opt_idx": <float>,
}
```

Where
- `"midpoint"` specifies the location of said boundary along the optical axis
  (`"Spherical"` specifies the position as `"midpoint" - "radius"`)
- `"height"` specifies the maximum absolute distance between the ray and an intersection
  point on the boundary that is possible
- `"opt_idx"` specifies the new optical index after the boundary

`Line` boundaries are defined by a point on the optical axis (`midpoint` / `m`) and a
maximum allowed distance (`height`) from the optical axis.
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=\color{orange}{x = m}">
</p>

`Spherical` boundaries are defined by a position (`midpoint` / `m`) on the optical axis
and a radius.
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=\color{orange}{(x%20-%20m)^2%20%2B%20y^2 = R^2}">
</p>

`Conic` boundaries are defined by the conic formula with all polynomial coefficients
set to zero.
<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=\color{orange}{z(r)%20=%20m%20%2B%20\frac{r^2}{R%20\left(1%20%2B%20\sqrt{1%20-%20(1%20%2B%20\kappa)%20\frac{r^2}{R^2}}\right)}}">
</p>

All of the boundary types have the additional restriction:
  <img src="https://render.githubusercontent.com/render/math?math=\color{orange}{y \leq \text{height}}">, constraining the distance of intersection points to the optical axis.
