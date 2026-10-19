# Validation record

Environment: Python 3.12, NumPy 2.3.5, SciPy 1.18.1, Matplotlib 3.11.2.
Baseline: NA=ND=1e16 cm^-3, P/N contact distances=2 um each,
tanh junction transition=20 nm, T=300 K, equal midgap SRH lifetimes=10 ns.

```bash
python run_all.py
python -m unittest discover -s tests -v
```

13 tests pass: stable Bernoulli identities; exact equilibrium zero flux;
contact neutrality and quasi-Fermi signs; independent collocation; continuity
sources and flux reconstruction; positive carriers; forward recombination and
reverse generation; no-SRH analytical limit; mesh order; bias path and mixing;
lifetime effect; asymmetric diode; temperature trend and input guards.

Main mesh: 601 nodes; main I-V: 25 biases. Three lifetime scans solve 31 biases
each, retaining 87 interior ideality-factor records. The inverse mesh has 301
nodes and the synthetic data contain 12 positive forward-bias observations.

| Quantity | Approximate result |
| --- | ---: |
| Max relative J difference at five collocation-check biases | 7.20e-5 |
| Max potential difference | 7.71e-6 V |
| Max normalized Poisson residual | 1.79e-10 |
| Max flux/activity log reconstruction discrepancy | 9.14e-10 |
| J at 0.20 V | 0.139062 A/m^2 |
| Fitted common lifetime, truth=8 ns | 7.99294 ns |
| Synthetic lifetime relative error | -0.0883% |

The reference uses continuous ODE collocation with the SG profile as its initial
guess. Its RMS residual is an equation residual, not a strict error bound.
Current constancy is structural in the integral SG evaluation and is not counted
as independent validation. Raw local currents remain available to diagnose
cancellation at low bias.

Files: validation.json, independent_reference.csv, mesh_convergence.csv,
gummel_traces.json, iv.csv, current_components.csv, lifetime_sweep.csv,
extraction.json and data/synthetic_lifetime_iv.csv (all under figures/).
The report describes this default baseline; custom cases require new records.
