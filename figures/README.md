# Reproducible transport results

Run `python run_all.py` from the repository root. Eight main PNG plots cover
I-V, carrier/quasi-Fermi profiles, bands, current continuity, lifetime/ideality,
independent collocation and mesh checks, solver diagnostics, and synthetic fitting.

CSV headers include units. J is a current density in A/m^2; multiply by device
area in m^2 to obtain A. R is a net rate in m^-3 s^-1, negative for generation.
Quasi-Fermi energies have the signs Fn=VT*a and Fp=-VT*b in eV.

`data/` contains explicitly synthetic, seeded observations and metadata. The
forward generator uses adaptive continuous collocation and the inverse uses
SG/Gummel. Other physical parameters are assumed known in the one-parameter fit.

Total-current flatness alone is structural, not an independent correctness
check. Raw local currents, activity reconstruction, mesh errors and collocation
diagnostics are included. See ../docs/MODEL.md for their meanings and limitations.
