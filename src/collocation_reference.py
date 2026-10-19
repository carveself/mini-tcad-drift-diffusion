"""Independent continuous ODE collocation check of the SG/Gummel solution.

solve_bvp solves six first-order equations with adaptive collocation, not the
SG discrete equations. An interpolated SG solution is an initial guess only.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_bvp

from dd_solver import Diode,Q,EPS,LREF


def reference(d:Diode,state:dict,tolerance:float=2e-6):
    voltage=state["voltage_V"]
    y=state["x"]/LREF
    # Scale current unknowns for conditioning; densities use ni, not their logs.
    # Log densities create artificial, extremely thin layers at ohmic contacts.
    jscale=max(abs(state["current_density_A_m2"]),.01)
    u=state["u"]
    initial=np.vstack((u,np.gradient(u,y,edge_order=2),
                       state["n"]/d.ni,state["p"]/d.ni,
                       np.interp(y,(y[:-1]+y[1:])/2,state["Jn"])/jscale,
                       np.interp(y,(y[:-1]+y[1:])/2,state["Jp"])/jscale))
    cu,ca,cb=d.contacts(voltage)
    left=np.array([cu[0],np.exp(cu[0]+ca[0]),np.exp(-cu[0]+cb[0])])
    right=np.array([cu[1],np.exp(cu[1]+ca[1]),np.exp(-cu[1]+cb[1])])
    def ode(pos,z):
        n=d.ni*z[2];p=d.ni*z[3]
        if d.recombination:
            rate=d.ni**2*(z[2]*z[3]-1)/(d.lifetime_s*(n+p+2*d.ni))
        else: rate=np.zeros_like(n)
        return np.vstack((z[1],-Q*LREF**2/(EPS*d.vt)*(p-n+d.doping(pos*LREF)),
                          z[2]*z[1]+LREF*jscale*z[4]/(Q*d.dn*d.ni),
                          -z[3]*z[1]-LREF*jscale*z[5]/(Q*d.dp*d.ni),
                          Q*LREF*rate/jscale,-Q*LREF*rate/jscale))
    def boundaries(za,zb):
        return np.r_[za[[0,2,3]]-left,zb[[0,2,3]]-right]
    result=solve_bvp(ode,boundaries,y,initial,tol=tolerance,max_nodes=15000)
    if not result.success:
        raise RuntimeError(f"Independent collocation failed: {result.message}")
    values=result.sol(y)
    current=float(np.mean((values[4]+values[5])*jscale))
    if np.any(values[2:4]<=0):
        raise RuntimeError("Collocation produced nonpositive carrier densities")
    return {"psi":values[0]*d.vt,"n":d.ni*values[2],"p":d.ni*values[3],
            "Jn":values[4]*jscale,"Jp":values[5]*jscale,
            "current_density_A_m2":current,"nodes":len(result.x),
            "max_rms_residual":float(np.max(result.rms_residuals)),
            "boundary_residual_inf":float(np.max(np.abs(boundaries(result.y[:,0],result.y[:,-1])))),
            "status":result.message}
