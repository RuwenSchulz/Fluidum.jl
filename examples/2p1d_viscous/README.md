# Fluidum 2+1D viscous hydrodynamics — worked examples

Three runnable setups for `matrix2d_visc!` (`src/Matrix/2d_viscous.jl`), the transverse-Cartesian,
boost-invariant Milne solver for

```
phi = (T, u^x, u^y, pi^{yy}, pi^{zz}, pi^{xy}, Pi)
```

Each produces a figure in `figures/`. They are written to be **copied and edited**.

```sh
julia --project=Julia Julia/Fluidum.jl/examples/2p1d_viscous/01_first_run.jl
julia --project=Julia Julia/Fluidum.jl/examples/2p1d_viscous/02_resolution_and_numerical_viscosity.jl
julia --project=Julia Julia/Fluidum.jl/examples/2p1d_viscous/03_elliptic_fireball.jl
```

| | what it shows |
|---|---|
| **01** the first run | the whole API, and the one configuration with a closed-form answer: a uniform state must integrate the 0+1D Israel–Stewart equations, which it does to 1e-13 |
| **02** resolution and numerical viscosity | the sound attenuation against the *exact* MIS dispersion relation, and the measurement that the first-order upwind scheme's own diffusivity is `D_num ≈ 0.75 Δx` — comparable to the physical `Γ_s` at production resolution |
| **03** an elliptic fireball | the sector ladder on a deformed Woods–Saxon, the causal envelope measured from the characteristic matrix, and where the run actually leaves it |

## Conventions, and the one that has caught people

`pi^{zz}` is the **mixed** Milne component `pi^eta_eta = tau^2 pi^{eta eta}`, not `pi^{eta eta}`;
`pi^{yy}` is the Cartesian `pi^y_y`. Same mixed convention as `viscous_1d`'s
`piphiphi`/`pietaeta`. The Navier–Stokes fixed points on a Bjorken background are

```
pi^eta_eta -> -4 eta/(3 tau)      pi^y_y -> +2 eta/(3 tau)      Pi -> -zeta/tau
```

and the relaxation is **plain Israel–Stewart** — no `delta_pipi`, no `tau_pipi`, no `lambda_piPi`.
Example 01 pins all of that to 1e-13; a `delta_pipi = (4/3) tau_pi` term of the kind FiVoHydro
carries would move `pi^eta_eta` by ~8 % on the state it reaches, and it is measurably absent. That
difference is why a Fluidum-vs-FiVoHydro comparison must be run with FiVo's `delta_pipi` switched
off to be like-for-like (`Julia/Projects/FiVoFluidumComparison/COMPARISON_2P1D.md` §3).

## Use `FLUIDUM_2D_DERIVED=1` for new work

The shipped generated matrix has ten wrong entries in its two shear rows: the `pi^{xy}` row is
driven by the one-sided velocity gradient `-2 eta d_x u^y` where Israel–Stewart asks for
`-2 eta sigma^{xy} = -eta (d_x u^y + d_y u^x)`. Consequences:

* it is **exact for irrotational transverse flow**, which is why smooth fireballs started from rest
  have looked right — including the published code-to-code comparison;
* it is wrong by `-2 eta omega_xy` wherever the flow carries vorticity: a rigid rotation, which must
  generate no shear at all, generates `pi^{xy}` at the full Navier–Stokes rate;
* it makes the characteristic structure **anisotropic** — the transverse shear channel propagates at
  `sqrt(2 C_s)` along x and at 0 along y.

`FLUIDUM_2D_DERIVED=1` (or `Fluidum.VISC_2D_DERIVED[] = true`) routes every 2+1D kernel, the
10-field `matrix2d_visc_HQ_BG!` included, to the matrix re-derived from scratch by
`Julia/tools/derive_2p1d_viscous.wls`. It is not slower. It is default **off** only so that results
already on disk keep their meaning. Full account: the note at the top of `src/Matrix/2d_viscous.jl`.

## Where the gates live

These are examples, not gates.

| | |
|---|---|
| `Julia/Fluidum.jl/test/test_2p1d_viscous.jl` | the fast testset, run by `Pkg.test()` — closed-form characteristic speeds, the isotropy assertions, the Bjorken limit |
| `Julia/Fluidum.jl/bench/bench_2p1d_viscous.jl` | the solver against closed forms: Bjorken, the exact sound dispersion, ideal Gubser, 2-D against the 1-D cylindrical solver |
| `Julia/Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl` | the derivation side: which matrix entries are wrong, and the corrected matrix's gates |
