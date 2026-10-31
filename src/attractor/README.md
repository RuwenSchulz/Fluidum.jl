# src/attractor

What this branch carries on top of upstream Fluidum. Every other file under `src/` and `test/` is
the upstream code plus the upstream pull request, apart from the repairs and the system listed below. Outside `src/` the branch also
has `deriv/` (the symbolic computations behind the generated code), `examples/`, `bench/`, and a wider
`Integrals` compat in `Project.toml` (phd-git pins an older SciML stack).

- `attractor.jl` holds the diagnostic switches (`HQ_TAUN_SCALE`, `FLUIDUM_HQ_MOMENTUM_DIMS`,
  `FLUIDUM_HQ_M2_MIXTURE`), the charm transport helpers phd-git calls (`hq_taun`, `hq_transport`,
  `hq_m2_mixture`), the `HQdiffusion` constructors, and the closure forms of the ideal charm matrices.
  It also accepts the earlier second-moment keys in the dict initializer.
- `profiles_cent1_cent2.jl` holds the `(cent1, cent2)` methods of `Profiles`.
- `detector_struct.jl` holds the detector table and this branch's `InitialParameters`.

`HQdiffusion(slope, offset, mass)` freezes the linear law below T_fo = 0.156 GeV, the LP1 transport
contract. A bare `LinearDiffusion(slope, offset, mass)` freezes it at upstream's 0.180.

## Upstream functions repaired on this branch (2026-10-30), not in the upstream pull request

- `matrix1d_visc`: `pressure(ϕ[1], Val(1), params.eos)` had no method.
- `one_d_ideal_matrix_ruwen6`: column 3 of rows 1-2 was exchanged between A_t and A_r.
- `matrix1d_ideal_HQ_br_p_fugacity!`: `n`, `dn_dmu`, `dn_dT` were never bound, and the derivatives of n went into the P_h slots.
- `matrix1d_ideal_HQ_br_p_density!`: the same P_h slots, and α inverted with the fm⁻³ ↔ GeV³ factor off. Both now use
  `charm_pressure` and `charm_alpha`.
- `thermodynamic(T, α, list)` of a charm list: the Hessian was incomplete (∂²_T n and ∂²_α n were 0).
- `matrix2d_visc_HQ!` / `two_d_viscous_HQ_matrix`: the ten wrong medium entries, `ϕ[6]` read as the fugacity, the entropy
  in the pressure slot, and κ passed where D_s is expected.
- `Profiles(x::TabulatedData; exp_tail = true)`: `exponential_tail_pointlike` needs `offset`; it is now a keyword, default 0.
- `hq_pressure`: used `x` and `ccbar`, which were not defined there. It is now the thermodynamics of `charm_quark_list`.

`two_d_viscous_matrix` keeps its ten wrong entries: `legacy = true` exists to reproduce it, and `deriv/derive_upstream.wls`
lists them.

## The fully coupled system

`matrix1d_two_way_cM!`, `matrix2d_two_way_cM!` (Matrix/1d_two_way_cM_cilindrical.jl, 2d_two_way_cM.jl), and
`two_way_problem(...; cM = true)`: the second moment in T^{μν}, with c_M = D_s/T, on the two-way layouts. Also not in the upstream
pull request.
