"""Steady 1D silicon drift-diffusion diode, SG fluxes and Gummel iteration.

Unknowns: u=psi/VT, a=ln(n/ni)-u, b=ln(p/ni)+u. The latter are
logarithms of quasi-Fermi activities. x increases from the P to N contact.
All physical quantities use SI; CLI dopings are cm^-3 and lengths um.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.linalg import solve_banded

Q = 1.602176634e-19
KB = 1.380649e-23
EPS = 11.7*8.8541878128e-12
LREF = 1e-6


@dataclass(frozen=True)
class Diode:
    na_cm3: float = 1e16
    nd_cm3: float = 1e16
    p_length_um: float = 2.
    n_length_um: float = 2.
    transition_nm: float = 20.
    temperature: float = 300.
    lifetime_s: float = 1e-8
    recombination: bool = True

    def __post_init__(self):
        vals = (self.na_cm3, self.nd_cm3, self.p_length_um, self.n_length_um,
                self.transition_nm, self.temperature, self.lifetime_s)
        if not all(np.isfinite(vals)):
            raise ValueError("All diode parameters must be finite")
        if not 1e15 <= min(self.na_cm3, self.nd_cm3) or max(self.na_cm3, self.nd_cm3)>1e17:
            raise ValueError("Use NA, ND in 1e15..1e17 cm^-3")
        if not .5 <= min(self.p_length_um, self.n_length_um) or max(self.p_length_um, self.n_length_um)>10:
            raise ValueError("Each neutral-base/contact length must be .5..10 um")
        if not 10 <= self.transition_nm <= 100 or not 250 <= self.temperature <= 400:
            raise ValueError("Use transition=10..100 nm and T=250..400 K")
        if not 1e-10 <= self.lifetime_s <= 1e-5:
            raise ValueError("Use lifetime=1e-10..1e-5 s")

    @property
    def na(self): return self.na_cm3*1e6
    @property
    def nd(self): return self.nd_cm3*1e6
    @property
    def vt(self): return KB*self.temperature/Q
    @property
    def ni(self):
        def eg(t): return 1.17-4.73e-4*t*t/(t+636.)
        return float(1e16*(self.temperature/300.)**1.5*np.exp(
            eg(300.)/(2*KB/Q*300.)-eg(self.temperature)/(2*self.vt)))
    @property
    def dn(self): return .135*self.vt
    @property
    def dp(self): return .048*self.vt
    @property
    def vbi(self): return self.vt*(np.arcsinh(self.na/(2*self.ni))+np.arcsinh(self.nd/(2*self.ni)))

    def doping(self, x):
        return (self.nd-self.na)/2+(self.nd+self.na)/2*np.tanh(np.asarray(x)/(self.transition_nm*1e-9))

    def contacts(self, voltage):
        up = -np.arcsinh(self.na/(2*self.ni))
        un = np.arcsinh(self.nd/(2*self.ni))-voltage/self.vt
        return np.array([up, un]), np.array([0., voltage/self.vt]), np.array([0., -voltage/self.vt])

    def summary(self):
        return {"NA_cm-3":self.na_cm3,"ND_cm-3":self.nd_cm3,
                "P_length_um":self.p_length_um,"N_length_um":self.n_length_um,
                "transition_nm":self.transition_nm,"temperature_K":self.temperature,
                "lifetime_s":self.lifetime_s,"SRH_enabled":self.recombination,
                "ni_cm-3":self.ni/1e6,"Vbi_V":float(self.vbi),
                "electron_mobility_cm2_Vs":1350.,"hole_mobility_cm2_Vs":480.}


def bernoulli(z):
    """B(z)=z/(exp(z)-1), stable at zero and large positive/negative z."""
    z = np.asarray(z, dtype=float)
    out = np.empty_like(z)
    small, pos, neg = np.abs(z)<1e-4, z>50, z< -50
    mid = ~(small|pos|neg)
    t = z[small]
    out[small] = 1-t/2+t*t/12-t**4/720+t**6/30240
    out[pos] = z[pos]*np.exp(-z[pos])/(1-np.exp(-z[pos]))
    out[neg] = -z[neg]/(1-np.exp(z[neg]))
    out[mid] = z[mid]/np.expm1(z[mid])
    return out


def grid(d: Diode, nodes: int):
    if not isinstance(nodes,int) or nodes<101 or nodes%2!=1:
        raise ValueError("Use an odd node count >=101")
    t = np.linspace(-1,1,nodes)
    s = np.sign(t)*np.expm1(2*np.abs(t))/np.expm1(2.)
    return s*np.where(t<0,d.p_length_um,d.n_length_um)*1e-6


def carriers(d,u,a,b):
    ln, lp = u+a, -u+b
    if max(np.max(np.abs(ln)),np.max(np.abs(lp)))>120:
        raise ValueError("Trial carrier exponent is outside safe classical range")
    return d.ni*np.exp(ln), d.ni*np.exp(lp)


def srh(d,u,a,b):
    if not d.recombination: return np.zeros_like(u)
    n,p = carriers(d,u,a,b)
    # np-ni^2 is evaluated through activities, avoiding cancellation at equilibrium.
    return d.ni**2*np.expm1(a+b)/(d.lifetime_s*(n+p+2*d.ni))


def poisson(d,x,u,a,b,voltage):
    y = x/LREF
    h = np.diff(y)
    vol = (h[:-1]+h[1:])/2
    norm = max(d.na,d.nd)
    alpha = Q*norm*LREF**2/(EPS*d.vt)
    dop = d.doping(x)/norm
    bc = d.contacts(voltage)[0]
    u = np.array(u,copy=True)
    u[[0,-1]] = bc
    def equation(v, jac=False):
        n,p = carriers(d,v,a,b)
        f = np.empty_like(v)
        f[0],f[-1] = v[0]-bc[0],v[-1]-bc[1]
        f[1:-1] = (v[1:-1]-v[:-2])/h[:-1]+(v[1:-1]-v[2:])/h[1:]-alpha*(p[1:-1]/norm-n[1:-1]/norm+dop[1:-1])*vol
        if not jac: return f
        band = np.zeros((3,len(v)))
        band[1,[0,-1]]=1
        band[1,1:-1] = 1/h[:-1]+1/h[1:]+alpha*(p[1:-1]+n[1:-1])/norm*vol
        band[0,2:] = -1/h[1:]
        band[2,:-2] = -1/h[:-1]
        return f,band
    for _ in range(80):
        f,band = equation(u,True)
        residual = float(np.max(np.abs(f)))
        if residual<2e-10: return u,residual
        delta=solve_banded((1,1),band,-f)
        step=min(1.,3./max(np.max(np.abs(delta)),1e-30))
        for _ in range(30):
            trial=u+step*delta
            if np.max(np.abs(trial))<100 and np.max(np.abs(equation(trial)))<(1-1e-4*step)*residual:
                u=trial
                break
            step/=2
        else: raise RuntimeError("Nonlinear Poisson line search failed")
    raise RuntimeError("Nonlinear Poisson did not converge")


def conductances(d,x,u):
    du=np.diff(u)
    h=np.diff(x)
    # Jn=g_n*(A_R-A_L), Jp=-g_p*(B_R-B_L), A=exp(a), B=exp(b).
    gn=Q*d.dn*d.ni*np.exp(u[:-1])*bernoulli(-du)/h
    gp=Q*d.dp*d.ni*np.exp(-u[:-1])*bernoulli(du)/h
    return gn,gp


def continuity_solve(g,sink,source,bc):
    count=len(sink)
    band=np.zeros((3,count))
    rhs=np.array(source,copy=True)
    diag=g[:-1]+g[1:]+sink[1:-1]
    band[1,[0,-1]]=1
    band[1,1:-1]=1
    band[0,2:]=-g[1:]/diag
    band[2,:-2]=-g[:-1]/diag
    rhs[1:-1]/=diag
    rhs[[0,-1]]=bc
    result=solve_banded((1,1),band,rhs)
    if np.any(result<=0) or not np.all(np.isfinite(result)):
        raise RuntimeError("Positive carrier activity lost in continuity solve")
    return result


def evaluate_currents(d,x,u,a,b,voltage):
    """Fluxes from the integrated discrete continuity/contact constraints.

    This is a stable evaluation of the same SG system, not a fitted diode law.
    It avoids subtracting nearly equal large activities in neutral regions.
    Raw edge differences are also returned as an independent diagnostic.
    """
    gn,gp=conductances(d,x,u)
    h=np.diff(x)
    volume=(h[:-1]+h[1:])/2
    rate=srh(d,u,a,b)
    cum=np.r_[0.,np.cumsum(Q*rate[1:-1]*volume)]
    rn,rp=1/gn,1/gp
    jn_left=(np.expm1(voltage/d.vt)-np.dot(cum,rn))/np.sum(rn)
    jp_left=(-np.expm1(-voltage/d.vt)+np.dot(cum,rp))/np.sum(rp)
    jn,jp=jn_left+cum,jp_left-cum
    raw_n=gn*np.exp(a[:-1])*np.expm1(np.diff(a))
    raw_p=-gp*np.exp(b[:-1])*np.expm1(np.diff(b))
    # Integrated contact drops check the reconstructed flux constitutive relation.
    return jn,jp,raw_n,raw_p,rate


def solve(d: Diode, voltage: float, nodes: int=601, initial: dict|None=None,
          mixing: float=.7, tolerance: float=2e-9, max_iterations: int=500):
    if not np.isfinite(voltage) or not -.3<=voltage<=.65:
        raise ValueError("Use -.3<=voltage<=.65 V; no breakdown/tunneling model")
    if not 0<mixing<=1: raise ValueError("mixing must be in (0,1]")
    x=grid(d,nodes)
    bc_u,bc_a,bc_b=d.contacts(voltage)
    if initial is None:
        a=np.linspace(*bc_a,nodes);b=np.linspace(*bc_b,nodes)
        u=np.arcsinh(d.doping(x)/(2*d.ni))-voltage/d.vt*(x-x[0])/(x[-1]-x[0])
    else:
        if len(initial["x"])!=nodes or not np.allclose(initial["x"],x,rtol=1e-12,atol=1e-20):
            raise ValueError("Continuation requires the same device grid")
        shift=(voltage-initial["voltage_V"])/d.vt*(x-x[0])/(x[-1]-x[0])
        u=initial["u"].copy()-shift
        a=initial["a"].copy()+shift
        b=initial["b"].copy()-shift
    volume=np.r_[0.,(np.diff(x)[:-1]+np.diff(x)[1:])/2,0.]
    trace=[]
    for iteration in range(max_iterations):
        old=np.vstack((u,a,b))
        u,pres=poisson(d,x,u,a,b,voltage)
        if abs(voltage)<1e-15:
            a[:]=0;b[:]=0
            break
        gn,gp=conductances(d,x,u)
        n,p=carriers(d,u,a,b)
        if d.recombination:
            den=d.lifetime_s*(n+p+2*d.ni)
            sink_n=Q*p/den*d.ni*np.exp(u)*volume
            source=Q*d.ni**2/den*volume
        else: sink_n=np.zeros(nodes);source=np.zeros(nodes)
        anew=np.log(continuity_solve(gn,sink_n,source,np.exp(bc_a)))
        nnew=d.ni*np.exp(u+anew)
        if d.recombination:
            den=d.lifetime_s*(nnew+p+2*d.ni)
            sink_p=Q*nnew/den*d.ni*np.exp(-u)*volume
            source=Q*d.ni**2/den*volume
        else: sink_p=np.zeros(nodes)
        bnew=np.log(continuity_solve(gp,sink_p,source,np.exp(bc_b)))
        a=(1-mixing)*a+mixing*anew
        b=(1-mixing)*b+mixing*bnew
        a[[0,-1]]=bc_a;b[[0,-1]]=bc_b
        change=float(np.max(np.abs(np.vstack((u,a,b))-old)))
        trace.append({"iteration":iteration,"max_update":change,"poisson_residual":pres})
        if change<tolerance: break
    else: raise RuntimeError(f"Gummel did not converge at V={voltage:g}; change={change:.3e}")
    u,pres=poisson(d,x,u,a,b,voltage)
    n,p=carriers(d,u,a,b)
    jn,jp,raw_n,raw_p,rate=evaluate_currents(d,x,u,a,b,voltage)
    current=float(np.mean(jn+jp))
    rawerr=max(float(np.max(np.abs(raw_n-jn))),float(np.max(np.abs(raw_p-jp))))
    reconstruction_error=0.
    for g,j,activity,bc,sign in zip(conductances(d,x,u),(jn,jp),(a,b),
                                   d.contacts(voltage)[1:],(1.,-1.)):
        increments=sign*j/g
        left=np.exp(bc[0])+np.r_[0.,np.cumsum(increments)]
        right=np.exp(bc[1])-np.r_[np.cumsum(increments[::-1])[::-1],0.]
        fraction=np.r_[0.,np.cumsum(1/g)]/np.sum(1/g)
        reconstructed=np.where(fraction<.5,left,right)
        if np.any(reconstructed<=0):
            raise RuntimeError("Flux-integrated activity lost positivity")
        reconstruction_error=max(reconstruction_error,float(np.max(np.abs(np.log(reconstructed)-activity))))
    return {"x":x,"u":u,"a":a,"b":b,"n":n,"p":p,"psi":u*d.vt,
            "Fn_eV":a*d.vt,"Fp_eV":-b*d.vt,
            "Jn":jn,"Jp":jp,"Jn_raw":raw_n,"Jp_raw":raw_p,"R":rate,
            "voltage_V":float(voltage),"current_density_A_m2":current,
            "iterations":iteration+1,"poisson_residual":pres,"trace":trace,
            "raw_constitutive_mismatch_A_m2":rawerr,
            "raw_mismatch_relative_current":rawerr/max(abs(current),1e-20),
            "flux_activity_reconstruction_log_error":reconstruction_error,
            "total_current_spread_A_m2":float(np.ptp(jn+jp))}


def sweep(d:Diode, voltages, nodes:int=601):
    state=solve(d,0.,nodes)
    results=[]
    for target in voltages:
        start=state["voltage_V"]
        steps=max(1,int(np.ceil(abs(target-start)/.025)))
        for v in np.linspace(start,target,steps+1)[1:]:
            state=solve(d,float(v),nodes,initial=state)
        results.append(state)
    return results
