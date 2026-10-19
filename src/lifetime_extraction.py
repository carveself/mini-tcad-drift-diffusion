"""Synthetic SRH lifetime extraction with independent collocation data."""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

from dd_solver import Diode,sweep
from collocation_reference import reference
from experiments import savefig,csvwrite,plt


def extract(out:Path,nodes:int=301):
    truth=Diode(lifetime_s=8e-9)
    voltage=np.linspace(.08,.3,12)
    initial=sweep(truth,voltage,401)
    clean=np.array([reference(truth,s)["current_density_A_m2"] for s in initial])
    seed=20261019;sigma=.02
    noisy=clean*np.exp(np.random.default_rng(seed).normal(0,sigma,len(clean)))
    evaluations=[]
    def residual(logtau):
        dev=replace(truth,lifetime_s=float(np.exp(logtau[0])))
        prediction=np.array([s["current_density_A_m2"] for s in sweep(dev,voltage,nodes)])
        if np.any(prediction<=0): raise RuntimeError("Positive forward current required for log fit")
        error=(np.log(prediction)-np.log(noisy))/sigma
        evaluations.append({"evaluation":len(evaluations),"lifetime_s":dev.lifetime_s,"weighted_sum_squared":float(error@error)})
        return error
    fit=least_squares(residual,np.log([2e-8]),bounds=(np.log([1e-10]),np.log([1e-6])),
                      diff_step=1e-4,ftol=1e-8,xtol=1e-8,gtol=1e-7,max_nfev=30)
    if not fit.success: raise RuntimeError(f"Lifetime fit failed: {fit.message}")
    tau=float(np.exp(fit.x[0]))
    prediction=np.array([s["current_density_A_m2"] for s in sweep(replace(truth,lifetime_s=tau),voltage,nodes)])
    metrics={"data_origin":"SYNTHETIC, not measured: continuous adaptive collocation forward data",
             "truth_lifetime_s":truth.lifetime_s,"fit_lifetime_s":tau,"relative_error":tau/truth.lifetime_s-1,
             "seed":seed,"log_noise_sigma":sigma,"inverse_mesh_nodes":nodes,
             "known_parameters":"Doping, geometry, temperature, mobility, ni and ideal ohmic contacts fixed; only common SRH lifetime fitted",
             "evaluations":len(evaluations),"fit_message":fit.message,
             "approx_relative_1sigma":float(np.sqrt((fit.fun@fit.fun/(len(voltage)-1))/(fit.jac.T@fit.jac)[0,0]))}
    data=out/"data";data.mkdir(parents=True,exist_ok=True)
    csvwrite(data/"synthetic_lifetime_iv.csv",[{"voltage_V":float(v),"clean_collocation_J_A_m2":float(c),
                "synthetic_noisy_J_A_m2":float(n),"SG_fitted_J_A_m2":float(p)} for v,c,n,p in zip(voltage,clean,noisy,prediction)])
    csvwrite(out/"fit_evaluations.csv",evaluations)
    (out/"extraction.json").write_text(json.dumps(metrics,indent=2)+"\n")
    (data/"metadata.json").write_text(json.dumps(metrics,indent=2)+"\n")
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    axes[0].semilogy(voltage,clean,label="Synthetic collocation truth")
    axes[0].semilogy(voltage,noisy,"o",ms=4,label="Synthetic lognormal-noise data")
    axes[0].semilogy(voltage,prediction,"--",label="Fitted SG/Gummel model")
    axes[0].set(xlabel="Forward voltage (V)",ylabel="Current density (A/m^2)",title="SRH lifetime extraction: synthetic data")
    axes[0].legend(fontsize=8)
    axes[1].plot(voltage,(np.log(prediction)-np.log(noisy))/sigma,"o-")
    axes[1].axhline(0,color="gray",ls=":")
    axes[1].set(xlabel="Forward voltage (V)",ylabel="Log residual / sigma",title=f"tau={tau:.4e} s; error={metrics['relative_error']:.2%}")
    for ax in axes: ax.grid(alpha=.25)
    savefig(fig,out,"lifetime_extraction.png")
    return metrics
