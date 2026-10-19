"""Small reproduction checkpoints for the transport project."""
import argparse
import json
from pathlib import Path
import numpy as np
from dd_solver import Diode,solve,sweep
from experiments import ROOT,generate
from lifetime_extraction import extract


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=("equilibrium","biased","experiments","extract"))
    parser.add_argument("--output",type=Path,default=ROOT/"reproduction")
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if args.stage in ("equilibrium","biased"):
        state=solve(Diode(),0.,301) if args.stage=="equilibrium" else sweep(Diode(),[.2],301)[0]
        np.savetxt(args.output/f"{args.stage}.csv",np.column_stack((state["x"],state["psi"],state["n"],state["p"],state["Fn_eV"],state["Fp_eV"])),
                   delimiter=",",header="x_m,psi_V,n_m-3,p_m-3,Fn_eV,Fp_eV",comments="")
        print(json.dumps({key:state[key] for key in ("voltage_V","current_density_A_m2","poisson_residual","flux_activity_reconstruction_log_error")},indent=2))
    elif args.stage=="experiments":print(json.dumps(generate(Diode(),args.output),indent=2))
    else:print(json.dumps(extract(args.output),indent=2))


if __name__=="__main__":main()
