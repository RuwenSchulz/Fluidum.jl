# Fluidum 1+1D viscous hydrodynamics — worked examples

Three runnable setups for `matrix1d_visc!` and `matrix1d_visc_HQ!` (`src/Matrix/1d_viscous_cilindrical.jl`,
`1d_viscous_HQ_cilindrical_fugacity.jl`), the cylindrical, boost-invariant Milne solvers for

```
phi = (T, u^r, pi^phi_phi, pi^eta_eta, Pi)              viscous_1d()
phi = (T, u^r, pi^phi_phi, pi^eta_eta, Pi, alpha, nu^r) HQ_viscous_1d()
```

Each produces a figure in `figures/`. They are written to be copied and edited.

The examples need `Plots` in the active environment. It is not a dependency of Fluidum, and nothing
in `src/`, `test/` or `bench/` needs it. From the package root:

```sh
julia --project=. examples/1p1d_viscous/01_first_run.jl
julia --project=. examples/1p1d_viscous/02_fireball_and_resolution.jl
julia --project=. examples/1p1d_viscous/03_charm_current.jl
```

| | what it shows | time |
|---|---|---|
| **01** the first run | the whole API (layout, discretisation, fluid properties, IC, driver, readout) on the one configuration with a closed-form answer: a radially uniform state must integrate the 0+1D Israel–Stewart equations | 19 s |
| **02** a fireball, the axis and resolution | the sector ladder (ideal → +shear → +bulk) from one initial condition; `u^r → 0` at the axis, faster under refinement, as the odd parity requires; and the self-convergence order, ≈1 because `basicupwinding` is first order | 26 s |
| **03** the charm current | `alpha` and `nu^r` on a cooling fireball, the current built from a gradient that starts at zero, and the charm-number integral over Cauchy slices, whose drift halves with each refinement | 100 s |

### 01 — the first run

A radially uniform state must integrate the 0+1D Israel–Stewart equations. The solver and the
closed form are drawn on top of each other.

![the Bjorken limit](figures/ex01_bjorken.png)

### 02 — a fireball, the axis, and resolution

Left: the sector ladder from one initial condition: ideal, +shear, +shear+bulk. Both dissipative
sectors resist, so both reduce the radial flow; a sector that increased it would be a sign error.
Right: `|u^r|` at the axis against resolution, which checks the odd parity.

![the fireball, the axis and resolution](figures/ex02_fireball_and_resolution.png)

### 03 — the charm current

`ν^r` starts at zero and is built by the gradient that cooling and radial flow create; it vanishes
at the axis by parity and where `α` is flat.

![the charm diffusion current](figures/ex03_charm_current.png)

## Conventions

`pietaeta` is the mixed Milne component `pi^eta_eta = tau^2 pi^{eta eta}`, not `pi^{eta eta}`;
`piphiphi` is the mixed `r^2 pi^{phi phi}`. `viscous_2d` uses the same mixed convention for
`piyy`/`pizz`. The Navier–Stokes fixed points on a Bjorken background are

```
pi^eta_eta -> -4 eta/(3 tau)      pi^phi_phi -> +2 eta/(3 tau)      Pi -> -zeta/tau
```

The relaxation is Israel–Stewart plus `delta_pipi pi theta`, with `delta_pipi = (4/3) tau_pi` (the
DNMR value) by default. `FLUIDUM_DELTA_PIPI=0` gives plain Israel–Stewart. There is no `tau_pipi`
and no `lambda_piPi`. Example 01's reference follows the switch; the two values differ by ~20 % in
`pi^eta_eta` on the state it reaches. When comparing with another code, match its `delta_pipi`.

## Pitfalls

1. **`r = 0` is a symmetry axis, not a boundary.** `viscous_1d` declares `ur` odd there and every
   other field even. A wrong parity does not crash; it puts a kink in the first few cells that
   looks like a physical feature. Example 02 measures `|u^r|` in the first cell against its peak and
   shows it falling monotonically under refinement; a parity error leaves it flat in `N`.

2. **The charm sector fails on a sharp fireball edge; the medium does not.** With the Woods–Saxon
   width `a = 0.6` of example 02, example 03 aborts at `tau ≈ 3.3–4.2` with *"dt was forced below
   floating point epsilon"*: `alpha` runs away (to about −729) where the charm density underflows,
   and the `alpha` row's time coefficient `u^tau dn_da = u^tau n` vanishes with it, leaving the
   `dmn_eps` regulator as the coefficient. Example 03 uses `a = 2.0` for that reason. The 2+1D
   solver fails the same way on realistic initial conditions, and alternative regulators do not
   help. The cause is the `alpha` runaway, not the regulator: a soft edge avoids it, it does not
   fix it.

3. **Julia soft scope.** `prev = drift` inside a top-level `for` creates a new local every
   iteration, so an accumulator silently never updates. Example 03 keeps its loop in a function.

## Reading the numbers

Quote `|pi|/e` above freeze-out only. In the dilute tail `e → 0` while `pi` relaxes on a finite
clock, so the ratio there is not a property of the fluid, and it runs past 1. Example 02 prints 0.08
with a `T > T_fo` mask and 1.96 without.

## Tests and benchmarks

| | |
|---|---|
| `test/test_1p1d_viscous.jl` | the fast testset, run by `Pkg.test()`: layout and parities, closed-form characteristic speeds, the shear closure, Bjorken, and the derived-matrix control |
| `bench/bench_1p1d_viscous.jl` | the solver against closed forms: the MIS front, the sound cone, `A_t` conditioning, the derivation control, Bjorken, the axis ladder and self-convergence (12 checks, ~25 s) |

```sh
julia --project=. bench/bench_1p1d_viscous.jl
```

The 2+1D examples are in [`../2p1d_viscous/`](../2p1d_viscous/README.md).
