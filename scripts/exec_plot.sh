#!/bin/bash

infile=${1:-infile.json}
cargo build
./target/debug/rs_raytrace $infile > path.csv
python ./scripts/plot_path.py $infile
