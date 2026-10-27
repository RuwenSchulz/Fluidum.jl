# src/attractor

What this branch carries on top of upstream Fluidum. Every other file under `src/` and `test/` is
the upstream code plus the open upstream pull request, byte for byte. Outside `src/` the branch also
has `deriv/` (symbolic derivations, for a later pull request), `examples/`, `bench/`, and a wider
`Integrals` compat in `Project.toml` (phd-git pins an older SciML stack).

- `attractor.jl` holds the diagnostic switches (`HQ_TAUN_SCALE`, `FLUIDUM_HQ_MOMENTUM_DIMS`,
  `FLUIDUM_HQ_M2_MIXTURE`), the charm transport helpers phd-git calls (`hq_taun`, `hq_transport`,
  `hq_m2_mixture`), the `HQdiffusion` constructors, and the closure forms of the ideal charm matrices.
  It also accepts the earlier second-moment keys in the dict initializer.
- `profiles_cent1_cent2.jl` holds the `(cent1, cent2)` methods of `Profiles`.
- `detector_struct.jl` holds the detector table and this branch's `InitialParameters`.

`HQdiffusion(slope, offset, mass)` freezes the linear law below T_fo = 0.156 GeV, the LP1 transport
contract. A bare `LinearDiffusion(slope, offset, mass)` freezes it at upstream's 0.180.

## Upstream functions found not to work

- `matrix1d_visc` (Matrix/1d_viscous_cilindrical.jl): `pressure(ϕ[1], Val(1), params.eos)` has no
  method, so it throws a MethodError on every EoS. `matrix1d_visc!` works.
- `matrix1d_ideal_HQ_br_p_fugacity!`: throws `UndefVarError: n`; `n`, `dn_dmu`, `dn_dT` are never bound.
- `matrix1d_ideal_HQ_br_p_density!`: builds a non-finite matrix at every state. The charm EoS returns
  an incomplete Hessian, so A_t is singular, and `invert_n_for_alpha` is off by the fm⁻³ ↔ GeV³ factor.
- `matrix2d_visc_HQ!` / `two_d_viscous_HQ_matrix` (Matrix/2d_viscous_fugacity.jl): complex and
  superluminal characteristic speeds on physical states; the wrapper reads `ϕ[6]` as the fugacity and
  passes the entropy into the pressure slot. Use `matrix2d_visc_HQ_BG!`.
- `Profiles(x::TabulatedData; exp_tail = true)` calls `exponential_tail_pointlike` without its
  required `offset` keyword.
- `hq_pressure` (EquationOfState/charm_HRG.jl) uses `x` and `ccbar`, which are not defined there.
