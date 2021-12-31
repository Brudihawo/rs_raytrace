use serde_json;

use std::fs::File;

mod util;
mod sim;
mod ray;
mod boundary;

use sim::{SimTemplate, RaySim};
use util::*;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        exit(
            format!(
                "Provide an infile name as the first argument!\n\
                 Usage:\n    \
                 {} <infile_name>",
                     args[0]
            ),
            1,
        );
    }

    let infile_name = &args[1];
    match File::open(infile_name) {
        Ok(infile) => {
            let res: Result<SimTemplate, serde_json::Error> = serde_json::from_reader(infile);
            match res {
                Ok(sim_temp) => {
                    let mut sim = RaySim::from_template(sim_temp);
                    sim.sim();
                }
                Err(error) => {
                    exit(
                        format!(
                            "Deserialisation of `{}` failed with error:\n{}!",
                            infile_name,
                            error.to_string()
                        ),
                        1,
                    );
                }
            }
        }
        Err(error) => {
            eprintln!(
                "Opening file {} failed with error: {}!",
                infile_name,
                error.to_string()
            );
            std::process::exit(1);
        }
    }
}
