#!/bin/bash

cargo build
./target/debug/rs_raytrace > path.csv
python ./scripts/plot_path.py
