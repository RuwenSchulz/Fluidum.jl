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

## A trap that bit this file's neighbour

`FLUIDUM_2D_DERIVED` was **dead on arrival** when first shipped, written as

```julia
const VISC_2D_DERIVED = Ref(get(ENV, "FLUIDUM_2D_DERIVED", "0") == "1")   # WRONG
```

at module top level. A `Ref` initialised at top level is evaluated when the **precompile image is
built** and its value serialised into the cache, so every later load carries the build-time value and
the environment variable does nothing. `get(ENV, ...)` inside the process returns `"1"` quite
happily; the `Ref` stays `false`. The fix is `Ref(false)` plus an assignment in `__init__`, which
runs on every module load — and that exact lesson was already written into `src/Fluidum.jl:58` for
`HQ_TAUN_SCALE`, three lines above its own fix, having been paid for once already.

⚠ **A dead switch and a genuinely null effect are indistinguishable from the output alone.** The A/B
that exposed it came back bit-identical in all 108 numbers, which reads exactly like "the repair
changes nothing". Before believing any null A/B on an env-controlled flag, assert the flag is
actually set inside the process:

```julia
julia -e 'using Fluidum; @show Fluidum.VISC_2D_DERIVED[]'                      # false
FLUIDUM_2D_DERIVED=1 julia -e 'using Fluidum; @show Fluidum.VISC_2D_DERIVED[]' # true
```

Gates in this file are immune by construction: M2 sets `F.VISC_2D_DERIVED[]` **programmatically**
and restores it, which always worked — only the `ENV` path was dead. `FLUIDUM_2DV_FULL` is read in a
*script*, which is not precompiled, so it is evaluated at run time and is fine.

## Companion

`Julia/Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl` carries the derivation side —
`Julia/tools/derive_2p1d_viscous.wls` regenerates the matrices from scratch, reproduces the shipped
1-D cylindrical production kernel entry by entry, and audits the shipped 2-D matrix against the
derived one. Read the two together: that file says which **entries** are wrong, this one says what
the **solver** does about it.
