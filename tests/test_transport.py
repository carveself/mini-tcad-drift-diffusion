"""Physical identities, independent collocation and numerical convergence."""
from pathlib import Path
import sys
import unittest
from dataclasses import replace
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from dd_solver import Diode,Q,bernoulli,conductances,solve,sweep
from collocation_reference import reference
from experiments import ideal_finite_base


class TransportTests(unittest.TestCase):
    def test_bernoulli_limits_and_symmetry(self):
        self.assertEqual(float(bernoulli(0.)),1.)
        x=np.array([-20.,-1.,-.00001,.00001,1.,20.])
        np.testing.assert_allclose(bernoulli(-x),np.exp(x)*bernoulli(x),rtol=1e-13)
        self.assertTrue(np.all(np.isfinite(bernoulli(np.array([-1000.,1000.])))))
        self.assertAlmostEqual(float(bernoulli(-1000.)),1000.)

    def test_equilibrium_exact_zero_flux(self):
        s=solve(Diode(),0.,301)
        self.assertEqual(s["current_density_A_m2"],0.)
        np.testing.assert_allclose(s["Jn_raw"],0,atol=0)
        np.testing.assert_allclose(s["Jp_raw"],0,atol=0)
        np.testing.assert_allclose(s["R"],0,atol=0)
        self.assertLess(s["poisson_residual"],2e-10)

    def test_contact_neutrality_and_fermi_sign(self):
        d=Diode();s=sweep(d,[.3],301)[0]
        np.testing.assert_allclose(s["n"][[0,-1]]-s["p"][[0,-1]],[-d.na,d.nd],rtol=1e-12)
        np.testing.assert_allclose(s["Fn_eV"][[0,-1]],[0,.3],atol=1e-14)
        np.testing.assert_allclose(s["Fp_eV"][[0,-1]],[0,.3],atol=1e-14)
        np.testing.assert_allclose(np.log(s["n"]*s["p"]/d.ni**2),(s["Fn_eV"]-s["Fp_eV"])/d.vt,atol=1e-12)

    def test_independent_collocation(self):
        d=Diode()
        for s in sweep(d,[-.2,.1,.4],301):
            r=reference(d,s)
            self.assertLess(abs(s["current_density_A_m2"]/r["current_density_A_m2"]-1),1e-3)
            self.assertLess(np.max(np.abs(s["psi"]-r["psi"])),1e-4)
            self.assertLess(r["max_rms_residual"],2.1e-6)

    def test_continuity_balance_with_srh(self):
        d=Diode();s=sweep(d,[.3],301)[0]
        h=np.diff(s["x"]);source=Q*s["R"][1:-1]*(h[:-1]+h[1:])/2
        np.testing.assert_allclose(np.diff(s["Jn"]),source,atol=1e-14,rtol=1e-8)
        np.testing.assert_allclose(np.diff(s["Jp"]),-source,atol=1e-14,rtol=1e-8)
        self.assertLess(s["total_current_spread_A_m2"],1e-12)
        self.assertLess(s["flux_activity_reconstruction_log_error"],5e-8)

    def test_positive_carriers(self):
        for s in sweep(Diode(),[-.2,.2,.55],301):
            self.assertTrue(np.all(s["n"]>0) and np.all(s["p"]>0))
            self.assertTrue(np.all(np.isfinite(s["psi"])))

    def test_srh_generation_and_recombination(self):
        d=Diode();rev,fwd=sweep(d,[-.2,.2],301)
        center=len(rev["x"])//2
        self.assertLess(rev["R"][center],0)
        self.assertGreater(fwd["R"][center],0)
        self.assertLess(rev["current_density_A_m2"],0)
        self.assertGreater(fwd["current_density_A_m2"],0)

    def test_no_recombination_low_injection_reference(self):
        d=Diode(recombination=False)
        for s in sweep(d,[.1,.2,.4],301):
            np.testing.assert_array_equal(s["R"],0)
            approx=float(ideal_finite_base(d,s["voltage_V"]))
            self.assertLess(abs(s["current_density_A_m2"]/approx-1),.02)

    def test_mesh_order_against_collocation(self):
        d=Diode();r=reference(d,sweep(d,[.45],301)[0])["current_density_A_m2"]
        err=[abs(sweep(d,[.45],n)[0]["current_density_A_m2"]/r-1) for n in (151,301,601)]
        self.assertTrue(3.4<err[0]/err[1]<4.4)
        self.assertTrue(3.4<err[1]/err[2]<4.4)

    def test_bias_path_and_mixing_independence(self):
        d=Diode();a=sweep(d,[-.2,.3],301)[-1];b=solve(d,.3,301,mixing=.35)
        self.assertLess(abs(a["current_density_A_m2"]/b["current_density_A_m2"]-1),1e-6)

    def test_lifetime_effect(self):
        short=sweep(Diode(lifetime_s=1e-9),[.2],301)[0]
        long=sweep(Diode(lifetime_s=1e-7),[.2],301)[0]
        self.assertGreater(short["current_density_A_m2"],long["current_density_A_m2"])

    def test_asymmetric_diode(self):
        d=Diode(na_cm3=3e16,nd_cm3=1e16)
        s=sweep(d,[.2],301)[0];r=reference(d,s)
        self.assertLess(abs(s["current_density_A_m2"]/r["current_density_A_m2"]-1),.002)

    def test_temperature_and_input_guards(self):
        self.assertLess(Diode(temperature=350).vbi,Diode(temperature=275).vbi)
        for kwargs in ({"na_cm3":0},{"lifetime_s":float("nan")},{"transition_nm":0}):
            with self.assertRaises(ValueError):Diode(**kwargs)
        with self.assertRaises(ValueError):solve(Diode(),.8)
        with self.assertRaises(ValueError):solve(Diode(),.2,nodes=200)
        with self.assertRaises(ValueError):solve(Diode(),.2,mixing=0)


if __name__=="__main__":unittest.main()
