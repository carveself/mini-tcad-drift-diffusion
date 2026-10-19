"""Transport experiments: SRH physics, mesh error, independent ODE checks."""
from dataclasses import replace
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from dd_solver import Diode,Q,EPS,sweep,solve
from collocation_reference import reference

ROOT=Path(__file__).resolve().parents[1]


def savefig(fig,out,name):
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/name,dpi=170,bbox_inches="tight")
    plt.close(fig)


def csvwrite(path,rows):
    with path.open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def ideal_finite_base(d,voltage):
    """Low-injection, no-SRH, depletion-edge finite-base approximation only."""
    voltage=np.asarray(voltage)
    width=np.sqrt(2*EPS/Q*(1/d.na+1/d.nd)*np.maximum(d.vbi-voltage,.01))
    xp=width*d.nd/(d.na+d.nd);xn=width*d.na/(d.na+d.nd)
    lp=d.p_length_um*1e-6-xp;ln=d.n_length_um*1e-6-xn
    if np.any(lp<=0) or np.any(ln<=0): raise ValueError("Neutral bases too short for analytical approximation")
    return Q*d.ni**2*(d.dn/(d.na*lp)+d.dp/(d.nd*ln))*np.expm1(voltage/d.vt)


def generate(d:Diode,out:Path,nodes:int=601):
    out.mkdir(parents=True,exist_ok=True)
    voltages=np.r_[[-.2,-.1,0.],np.linspace(.025,.55,22)]
    states=sweep(d,voltages,nodes)
    idealdev=replace(d,recombination=False)
    idealstates=sweep(idealdev,voltages,nodes)
    j=np.array([s["current_density_A_m2"] for s in states])
    j0=np.array([s["current_density_A_m2"] for s in idealstates])
    analytic=ideal_finite_base(d,voltages)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    for ax in axes:
        ax.plot(voltages,j,label="SG + SRH")
        ax.plot(voltages,j0,"--",label="SG, recombination disabled")
        ax.plot(voltages,analytic,":",label="Finite-base low-injection approximation (no SRH)")
        ax.set(xlabel="Applied forward voltage (V)",ylabel="Current density (A/m^2)")
        ax.grid(alpha=.25);ax.legend(fontsize=8)
    axes[0].set_title("Self-consistent drift-diffusion I-V")
    axes[1].set_yscale("symlog",linthresh=1e-5)
    axes[1].set_title("Reverse generation and forward injection")
    savefig(fig,out,"iv_comparison.png")
    rows=[]
    for s,without,approx in zip(states,idealstates,analytic):
        rows.append({"voltage_V":s["voltage_V"],"J_SRH_A_m2":s["current_density_A_m2"],
                     "J_no_R_A_m2":without["current_density_A_m2"],"J_analytic_approx_A_m2":float(approx),
                     "iterations_last_step":s["iterations"],"poisson_residual":s["poisson_residual"],
                     "current_spread_A_m2":s["total_current_spread_A_m2"],
                     "flux_reconstruction_log_error":s["flux_activity_reconstruction_log_error"],
                     "raw_SG_mismatch_relative_J":s["raw_mismatch_relative_current"]})
    csvwrite(out/"iv.csv",rows)

    profile_states=sweep(d,[-.2,0.,.2,.5],nodes)
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for s in profile_states:
        x=s["x"]*1e6;label=f"V={s['voltage_V']:.2f} V"
        axes[0,0].plot(x,s["psi"],label=label)
        axes[0,1].semilogy(x,s["n"]/1e6,label=label+", n")
        axes[0,1].semilogy(x,s["p"]/1e6,"--",alpha=.65)
        axes[1,0].plot(x,s["Fn_eV"],label=label+", Fn")
        axes[1,0].plot(x,s["Fp_eV"],"--",alpha=.65)
        axes[1,1].plot(x,s["R"],label=label)
        name=f"profile_{s['voltage_V']:+.2f}V.csv".replace("+","plus").replace("-","minus")
        np.savetxt(out/name,np.column_stack((s["x"],s["psi"],s["n"],s["p"],s["Fn_eV"],s["Fp_eV"],s["R"])),
                   delimiter=",",header="x_m,psi_V,n_m-3,p_m-3,Fn_eV,Fp_eV,R_m-3_s-1",comments="")
    titles=("Electrostatic potential","Carriers (solid n, dashed p)","Quasi-Fermi energies (solid Fn, dashed Fp)","SRH rate (negative means generation)")
    labels=("Potential (V)","Carrier density (cm^-3)","Energy relative to P contact (eV)","SRH net rate (m^-3 s^-1)")
    for ax,title,label in zip(axes.flat,titles,labels):
        ax.set(title=title,ylabel=label,xlabel="Position (um), P to N")
        ax.grid(alpha=.25);ax.legend(fontsize=7)
    savefig(fig,out,"profiles_quasi_fermi.png")
    fig,axes=plt.subplots(1,3,figsize=(14,4.7),constrained_layout=True)
    eg=1.17-4.73e-4*d.temperature**2/(d.temperature+636.)
    for ax,s in zip(axes,profile_states[1:]):
        x=s["x"]*1e6;ei=-s["psi"]
        ax.plot(x,ei+eg/2,label="Ec")
        ax.plot(x,ei-eg/2,label="Ev")
        ax.plot(x,s["Fn_eV"],"--",label="Fn: electron quasi-Fermi energy")
        ax.plot(x,s["Fp_eV"],":",label="Fp: hole quasi-Fermi energy")
        ax.set(title=f"Bias={s['voltage_V']:.2f} V",xlabel="Position (um)",ylabel="Energy relative to P-contact EF (eV)",xlim=(-.8,.8))
        ax.legend(fontsize=7);ax.grid(alpha=.25)
    savefig(fig,out,"bands_quasi_fermi.png")

    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    s=profile_states[2];xmid=(s["x"][:-1]+s["x"][1:])/2
    axes[0].plot(xmid*1e6,s["Jn"],label="Electron current")
    axes[0].plot(xmid*1e6,s["Jp"],label="Hole current")
    axes[0].plot(xmid*1e6,s["Jn"]+s["Jp"],"--",label="Total current")
    axes[0].set(title="Continuity-integrated SG fluxes, V=0.20 V",xlabel="Position (um)",ylabel="Current density (A/m^2)")
    axes[1].plot(xmid*1e6,s["Jn_raw"]+s["Jp_raw"],label="Direct local SG differences")
    axes[1].axhline(s["current_density_A_m2"],ls="--",color="black",label="Integral evaluation of same SG system")
    axes[1].set(title="Floating-point cancellation in neutral regions",xlabel="Position (um)",ylabel="Total current density (A/m^2)")
    for ax in axes: ax.legend(fontsize=8);ax.grid(alpha=.25)
    savefig(fig,out,"current_continuity.png")
    csvwrite(out/"current_components.csv",[{"x_mid_m":float(x),"Jn_A_m2":float(n),"Jp_A_m2":float(p),
              "Jn_raw_A_m2":float(rn),"Jp_raw_A_m2":float(rp)} for x,n,p,rn,rp in zip(xmid,s["Jn"],s["Jp"],s["Jn_raw"],s["Jp_raw"])])

    forward=np.linspace(.05,.55,31)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    lifetime_rows=[]
    for lifetime in (1e-9,1e-8,1e-7):
        dev=replace(d,lifetime_s=lifetime,recombination=True)
        results=sweep(dev,forward,nodes)
        current=np.array([s["current_density_A_m2"] for s in results])
        # Numerical derivative, endpoints omitted in display and data.
        factor=1/(dev.vt*np.gradient(np.log(current),forward,edge_order=2))
        axes[0].semilogy(forward,current,label=f"tau={lifetime:g} s")
        axes[1].plot(forward[1:-1],factor[1:-1],label=f"tau={lifetime:g} s")
        for v,cur,f in zip(forward[1:-1],current[1:-1],factor[1:-1]):
            lifetime_rows.append({"lifetime_s":lifetime,"voltage_V":float(v),"J_A_m2":float(cur),"apparent_ideality_factor":float(f)})
    axes[0].set(xlabel="Forward voltage (V)",ylabel="Current density (A/m^2)",title="SRH lifetime experiment")
    axes[1].axhline(1,ls=":",color="gray");axes[1].axhline(2,ls=":",color="gray")
    axes[1].set(xlabel="Forward voltage (V)",ylabel="n_app = 1 / (VT d ln J / dV)",title="Apparent ideality is bias-dependent")
    for ax in axes: ax.legend(fontsize=8);ax.grid(alpha=.25)
    savefig(fig,out,"lifetime_ideality.png")
    csvwrite(out/"lifetime_sweep.csv",lifetime_rows)

    reference_voltages=(-.2,.1,.2,.4,.55)
    checked=sweep(d,reference_voltages,nodes)
    comparison=[]
    for s in checked:
        ref=reference(d,s)
        comparison.append({"voltage_V":s["voltage_V"],"SG_J_A_m2":s["current_density_A_m2"],
                           "collocation_J_A_m2":ref["current_density_A_m2"],
                           "relative_current_error":s["current_density_A_m2"]/ref["current_density_A_m2"]-1,
                           "max_potential_error_V":float(np.max(np.abs(s["psi"]-ref["psi"]))),
                           "collocation_nodes":ref["nodes"],"collocation_RMS_residual":ref["max_rms_residual"]})
    csvwrite(out/"independent_reference.csv",comparison)
    refinements=[]
    refcurrent=reference(d,sweep(d,[.45],301)[0])["current_density_A_m2"]
    for count in (151,301,601,1201):
        state=sweep(d,[.45],count)[0]
        refinements.append({"nodes":count,"current_density_A_m2":state["current_density_A_m2"],
                            "relative_error_vs_collocation":abs(state["current_density_A_m2"]/refcurrent-1)})
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    axes[0].plot([r["voltage_V"] for r in comparison],np.abs([r["relative_current_error"] for r in comparison])*100,"o-")
    axes[0].set(xlabel="Applied voltage (V)",ylabel="Current difference (%)",title="Independent adaptive ODE collocation")
    nn=np.array([r["nodes"]-1 for r in refinements]);ee=np.array([r["relative_error_vs_collocation"] for r in refinements])
    axes[1].loglog(nn,ee,"o-",label="Current error vs collocation")
    axes[1].loglog(nn,ee[0]*(nn[0]/nn)**2,"--",label="Second-order guide")
    axes[1].set(xlabel="Mesh intervals",ylabel="Relative current error",title="Mesh study at V=0.45 V")
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.25)
    savefig(fig,out,"reference_mesh.png")
    csvwrite(out/"mesh_convergence.csv",refinements)

    direct=solve(d,.4,nodes)
    conservative=solve(d,.4,nodes,mixing=.35)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    for state,label in ((direct,"Gummel mixing=0.70"),(conservative,"Gummel mixing=0.35")):
        axes[0].semilogy([r["iteration"] for r in state["trace"]],[r["max_update"] for r in state["trace"]],label=label)
    axes[0].set(xlabel="Gummel outer iteration",ylabel="Max dimensionless update",title="From neutral-potential initial guess")
    axes[0].legend(fontsize=8)
    axes[1].semilogy(voltages,np.maximum([s["flux_activity_reconstruction_log_error"] for s in states],1e-16),"o-",label="Flux-to-activity reconstruction error")
    axes[1].semilogy(voltages,np.maximum([s["poisson_residual"] for s in states],1e-16),"o-",label="Poisson residual")
    axes[1].set(xlabel="Applied voltage (V)",ylabel="Dimensionless diagnostic",title="Constitutive consistency and Poisson")
    axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.25)
    savefig(fig,out,"solver_diagnostics.png")
    (out/"gummel_traces.json").write_text(json.dumps({"mix_070":direct["trace"],"mix_035":conservative["trace"]},indent=2)+"\n")
    metrics={"device":d.summary(),"nodes":nodes,"IV_points":len(states),
             "lifetime_solved_biases":3*len(forward),"lifetime_retained_ideality_points":len(lifetime_rows),
             "max_reference_relative_current_error":max(abs(r["relative_current_error"]) for r in comparison),
             "max_reference_potential_error_V":max(r["max_potential_error_V"] for r in comparison),
             "max_poisson_residual":max(s["poisson_residual"] for s in states),
             "max_flux_activity_reconstruction_log_error":max(s["flux_activity_reconstruction_log_error"] for s in states),
             "max_total_current_spread_A_m2":max(s["total_current_spread_A_m2"] for s in states),
             "raw_current_note":"Raw local differences suffer cancellation at low currents. Total current is obtained from the integrated SG/contact constraints, not a diode formula.",
             "conservation_note":"Exact total-flux constancy is structural, not independent proof of correctness; activity reconstruction, collocation and mesh checks are also required.",
             "independent_reference":comparison,"mesh":refinements,
             "mixing_current_relative_difference":abs(direct["current_density_A_m2"]/conservative["current_density_A_m2"]-1)}
    (out/"validation.json").write_text(json.dumps(metrics,indent=2)+"\n")
    return metrics
