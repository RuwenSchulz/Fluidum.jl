# Fluidum benchmarks

| | |
|---|---|
| `bench_2p1d_viscous.jl` | the 2+1D viscous solver (`matrix2d_visc!`) against closed forms |

```sh
julia --project=Julia Julia/Fluidum.jl/bench/bench_2p1d_viscous.jl
FLUIDUM_2DV_FULL=1 julia --project=Julia Julia/Fluidum.jl/bench/bench_2p1d_viscous.jl   # finer grids
```

20 gates, ~6 minutes in fast mode. Every reference is a closed form or an ODE integrated inside the
file; nothing reads the DataStore, so the suite is portable and deterministic.

| section | what it establishes |
|---|---|
| **M** | the shear closure is traceless and u-orthogonal; the Israel–Stewart front at rest is `sqrt(c_s^2 + (4/3) eta/(w tau_pi) + zeta/(w tau_Pi))` to 8e-16; the transverse shear channel is `sqrt(2 C_s)` independent of T and of the EoS; with viscosity off the extreme characteristics are the relativistic sound cone; the causal envelope, measured; `A_t` conditioning down to T = 0.02 GeV |
| **B** | a transversely uniform state integrates the 0+1D Israel–Stewart equations to 1e-13 — the gate that pins the conventions and the absence of `delta_pipi` |
| **S** | the sound phase velocity and attenuation against the **exact** MIS dispersion root, and the measurement that the excess damping is first order in `Δx`, i.e. the upwind scheme's own numerical viscosity, `D_num ≈ 0.75 Δx` |
| **G** | ideal Gubser evolved in 2-D, first-order convergence, x↔y symmetric to round-off |
| **A** | the 1-D and 2-D transport-coefficient paths agree; an axisymmetric 2-D viscous solve converges onto the 1-D cylindrical solver in all five fields; axisymmetry is preserved in the limit |
| **C** | cost |

## What this suite found

Gate **M2** — rotational isotropy of the characteristic speeds — fails on the shipped matrix by
**0.69 c**. The `pi^{xy}` row is driven by the one-sided velocity gradient rather than its symmetric
part, so the transverse shear channel propagates at `sqrt(2 C_s)` along x and at 0 along y. The
corrected matrix (`FLUIDUM_2D_DERIVED=1`) passes at 2.8e-14; M2 runs both and reports both.

Until 2026-09-03 this solver's entire test coverage was one `@test all(isfinite.(result))`.

## Companion

`Julia/Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl` carries the derivation side —
`Julia/tools/derive_2p1d_viscous.wls` regenerates the matrices from scratch, reproduces the shipped
1-D cylindrical production kernel entry by entry, and audits the shipped 2-D matrix against the
derived one. Read the two together: that file says which **entries** are wrong, this one says what
the **solver** does about it.
