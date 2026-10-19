"""Run 1D SG drift-diffusion experiments and synthetic lifetime extraction."""
from pathlib import Path
import argparse
import json
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"src"))
from dd_solver import Diode
from experiments import generate
from lifetime_extraction import extract


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes",type=int,default=601,help="Odd main mesh node count >=101")
    parser.add_argument("--fit-nodes",type=int,default=301)
    parser.add_argument("--output",type=Path,default=ROOT/"figures")
    args=parser.parse_args()
    try:
        metrics=generate(Diode(),args.output,args.nodes)
        fitted=extract(args.output,args.fit_nodes)
    except ValueError as exc: parser.error(str(exc))
    print(json.dumps({"max_reference_current_error":metrics["max_reference_relative_current_error"],
          "max_reference_potential_error_V":metrics["max_reference_potential_error_V"],
          "max_flux_activity_reconstruction_error":metrics["max_flux_activity_reconstruction_log_error"],
          "fit_lifetime_s":fitted["fit_lifetime_s"],"fit_relative_error":fitted["relative_error"]},indent=2))
    print(f"Complete. Figures and data saved to {args.output.resolve()}")


if __name__=="__main__":main()
