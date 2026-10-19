# two_way.jl: the conformal set and the coupled two-way medium + charm system.
#
# Wrappers around the generated matrices of two_way_generated.jl (deriv/derive_two_way.wls).
# Five systems:
#
#   S1 ideal hydro           S2 viscous hydro (IS, δ_ππ = dpp_factor·τ_π, default 4/3)
#   S3 charm first moment on a prescribed background (consistent closure)
#   S4 S3 + the consistent charm second moment on a prescribed background
#   S5 the two-way system: ∇_μ(T_m + T_Q)^{μν} = 0, medium π/Π rows and the charm rows of S4,
#      evolved together; c = 1 ⇔ c_M = D_s/T (δM in T_Q and in the current row).
#
#   0+1D Bjorken: an ODE right-hand side (`tw_bjorken_rhs`) and a driver (`two_way_bjorken`).
#   1+1D cylindrical / 2+1D Cartesian: closures in Fluidum's convention
#   (A_i, Source, ϕ, tau, X, params) that write A_t⁻¹A_i and A_t⁻¹src without allocating,
#   for `oneshoot` / `problem` (fluidevo/discreteFields.jl).
#
# Conventions:
#   * δ_ππ = 4/3 τ_π (dpp_factor), δ_ΠΠ = 0 (dPP_factor).
#   * The corotation (vorticity) coupling wv of the 2+1D charm second moment is on by default.
#   * Medium in the Landau frame; the charm energy is slaved to T.
#   * h in T_Q is the thermodynamic enthalpy per particle of the charm EOS, and the charm
#     relaxation time is τ_n = D_s h/T. HQ_TAUN_SCALE and FLUIDUM_HQ_MOMENTUM_DIMS act on
#     `hq_transport` only, not here. At their defaults both give the same τ_n.
#   * κ uses the unregulated dn/dα. The vacuum regulator below enters the number row's
#     α-column only.
#   * 7-field 2+1D charm second moment: π_Q^η_η is eliminated algebraically.
#
# The charm EOS is the keyword `charm`: a hadron list (`Heavy_Quark` or its
# `HadronResonaceGas_ccbar` hadron_list), or `BoltzmannCharm(m, g)` (one massive Boltzmann
# quark in closed form; `charm_quark_list` builds the same as a one-particle list).
#   n, ∂_Tn, ∂_αn       the list's `thermodynamic(T, α, list)`
#   h                   Σ q² n_i h_i / Σ q² n_i,  h_i = m_i K3/K2. With this weighting D_s h/T
#                       equals `τ_diffusion_hadron` (Eq. 30 of 2205.07692) for any list.
#   ∂_Th                analytic, at fixed α
#   I31 = n T h, bQ = Th. For one species aQ = m_q² + 5Th and cQ = h K4K2/K3² = m K4/K3
#     (z = m_q/T, m_q = diffusion.mass), so τ_M = D_s cQ/(2T) = τ_n K4K2/(2K3²) is the τ_M of
#     `hq_transport`. For a multi-species list cQ = Σq²I42_i C_i/Σq²I42_i and aQ = ⟨m²⟩ + 5Th
#     (`hq_m2_mixture`); `mixture = false` / FLUIDUM_HQ_M2_MIXTURE=0 uses the one-species forms.
#   The derivation assumes unit-charge Boltzmann species (P0 = nT, h = h(T)). For a doubly charmed
#   state (J/ψ, q = 2) P0 = nT over-counts its pressure by n_J/ψ T and ∂_αh ≠ 0 is dropped. Both
#   are O(e^{2α} n_J/ψ).
#
# The medium EOS is `params.eos`. `ConformalEOS` selects the conformal generated functions (no bulk
# field; aSB = `conformal_aSB(eos)`); any other EOS selects the generic ones. `generic = true`
# forces the generic path for a ConformalEOS (the two paths agree). Transport: `params.shear`,
# `params.bulk`, and D_s from `params.diffusion` (ConstDiffusion / LinearDiffusion: D_s = DsT/T).
#
# Vacuum tail. At n → 0 (α → -∞) the number row's A_t[α, α] = u^τ ∂_αn → 0 and A_t is singular.
#   (i)  `dmn_eps` (default 1e-8 fm⁻³) is added to ∂_αn in the number row's α-column only, so the
#        row conserves n + dmn_eps·α instead of n. Negligible where n ≫ dmn_eps.
#   (ii) α is clamped to ±α_max (default 200) in the thermodynamics only.
# Cold edge. Below the freeze-out temperature the charm hydro is outside its validity
# (Knudsen ≫ 1: |ν| > n, α → -10³), and without treatment charm that crosses T_fo is destroyed.
# The region T > T_fo does not feel it. Two parts, both on by default (`edge = true`):
#   (v)  Decoupling (`tw_edge_factor`, in every closure here): D_s → f(T) D_s, f = 1 for T ≥ T_edge
#        (0.156 GeV), a C² smoothstep down to `edge_floor` (0.01) at T_edge − edge_width (20 MeV).
#        κ, τ_n, τ_M and c_M = c D_s/T scale with it, so ν and δM relax to ~0 and frozen-out charm
#        comoves with the medium. The equations for T ≥ T_edge are unchanged.
#   (vi) Flux-form corona number update (`two_way_problem` only, since it needs the neighbours): see
#        the comment above `TwoWayCoronaRHS`. At an under-resolved charm edge the quasi-linear α
#        update is not positivity-preserving and decoupling alone is not enough.
#   Use `two_way_problem(d, params, init, tspan; charm, …)` for runs. The closures two_way_1p1d! /
#   two_way_2p1d! carry (v) only.
# Two optional extras, off by default and inert in the physical region:
#   (iii) `nu_cap` > 0: damp the current where its drift |ν|/n exceeds nu_cap (rate up to
#        (1 + nu_cap_rate)/τ_n). Unstable under Tsit5.
#   (iv) `cM_nfloor` > 0: c → c·n/(n + cM_nfloor), as HQC_CM_NFLOOR.

# ── the charm EOS ────────────────────────────────────────────────────────────────────────────

"""
    BoltzmannCharm(; m = 1.5, g = 6.0)

One massive Boltzmann charm quark (unit charge, degeneracy `g`, mass `m` in GeV) in closed form.
n = g m² T K₂(m/T) e^α/(2π²) (fm⁻³), h = m K₃/K₂. Identical to the one-particle list
`charm_quark_list`.
"""
struct BoltzmannCharm{T} <: EquationOfState
    m::T
    g::T
end
BoltzmannCharm(; m = 1.5, g = 6.0) = BoltzmannCharm(promote(float(m), float(g))...)

"""
    charm_quark_list(; m = 1.5, g = 6) -> HadronResonaceGas_ccbar

The single charm quark as a one-particle `HadronResonaceGas_ccbar` list, for runs that want the
list code path.
"""
function charm_quark_list(; m = 1.5, g = 6)
    p = Particle("Charm", float(m), 0.0, Int(g), 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4)
    return HadronResonaceGas_ccbar(StructArray([p]), 1.0)
end

_tw_charm_list(x::Heavy_Quark) = x.hadron_list
_tw_charm_list(x::HadronResonaceGas_ccbar) = x
_tw_charm_list(x::BoltzmannCharm) = x
_tw_charm_list(::Nothing) = nothing

# d/dz (K3/K2) from exponentially SCALED Bessels (the e^z cancels in the ratio)
@inline function _tw_dK32(K1, K2, K3, K4)
    dK2 = -(K1 + K3) / 2
    dK3 = -(K2 + K4) / 2
    return (dK3 * K2 - K3 * dK2) / K2^2
end

# n, ∂_Tn, ∂_αn, h, ∂_Th for one Boltzmann species of charge q (all scaled-Bessel, no overflow)
@inline function _tw_species(T, α, m, g, q)
    z = m / T
    K1 = besselk1x(z); K2 = besselkx(2, z)
    K3 = K1 + 4 / z * K2
    K4 = K2 + 6 / z * K3
    ex = exp(q * α - z)
    n0 = g * m^2 * T * K2 * ex / (2π^2) * fmGeV3                       # particle density
    dn0 = g * ex * m^2 * (m * K1 + 2 * T * K2 + m * K3) / (4π^2 * T) * fmGeV3
    h = m * K3 / K2
    dh = -(z / T) * m * _tw_dK32(K1, K2, K3, K4)
    return n0, dn0, h, dh
end

function _tw_charm_core(T, α, x::BoltzmannCharm)
    n, dn, h, dh = _tw_species(T, α, x.m, x.g, 1)
    return (n = n, dndT = dn, dnda = n, h = h, dhdT = dh)
end

function _tw_charm_core(T, α, x::HadronResonaceGas_ccbar)
    th = thermodynamic(T, α, x)
    n = th.pressure
    dndT, dnda = th.pressure_derivative
    W = zero(T); D = zero(T); dW = zero(T); dD = zero(T)
    for i in eachindex(x)
        q = round(x[i].Nc + x[i].Nac)
        m = x[i].Mass
        m / T < 1000 || continue
        n0, dn0, hi, dhi = _tw_species(T, α, m, float(x[i].Degeneracy), q)
        w = q^2 * n0; dw = q^2 * dn0
        W += w * hi; D += w
        dW += dw * hi + w * dhi; dD += dw
    end
    if D > 0
        h = W / D
        dhdT = (dW - h * dD) / D
    else                                    # every species underflowed: the lightest's ratio
        _, _, h, dhdT = _tw_species(T, zero(α), x[1].Mass, 1.0, 1)
    end
    return (n = n, dndT = dndT, dnda = dnda, h = h, dhdT = dhdT)
end

"""
    two_way_charm_thermo(T, α, charm, mq) -> NamedTuple

Every charm-EOS input of the generated functions at (T, α):
`n, dndT, dnda, h, dhdT, I31, dI31dT, dI31da, aQ, bQ, dbQdT, cQ, dcQdT` (header of
two_way.jl for the definitions; `mq` is the mass in aQ and in the K-ratio of cQ,
`params.diffusion.mass` in the wrappers). `charm` is a `BoltzmannCharm`, a
`HadronResonaceGas_ccbar` list or a `Heavy_Quark`.
"""
function two_way_charm_thermo(T, α, charm, mq; mixture::Bool = HQ_M2_MIXTURE[])
    L = _tw_charm_list(charm)
    c = _tw_charm_core(T, α, L)
    n, dndT, dnda, h, dhdT = c.n, c.dndT, c.dnda, c.h, c.dhdT
    if mixture && _tw_is_mixture(L, mq)
        # multi-species list: C, ∂_T C and the m² of A from the species sums (hq_m2_mixture,
        # shared with the m2 drivers)
        x = hq_m2_mixture(T, α, L, mq)
        return (n = n, dndT = dndT, dnda = dnda, h = h, dhdT = dhdT,
                I31 = n * T * h, dI31dT = dndT * T * h + n * h + n * T * dhdT, dI31da = dnda * T * h,
                aQ = x.m2bar + 5 * T * h, bQ = T * h, dbQdT = h + T * dhdT, cQ = x.C, dcQdT = x.dCdT)
    end
    z = mq / T
    K1 = besselk1x(z); K2 = besselkx(2, z)
    K3 = K1 + 4 / z * K2
    K4 = K2 + 6 / z * K3
    K5 = K3 + 8 / z * K4
    R = K4 * K2 / K3^2
    dlnR = -(K3 + K5) / (2 * K4) - (K1 + K3) / (2 * K2) + (K2 + K4) / K3     # d ln R / dz
    cQ = h * R
    dcQdT = dhdT * R + h * R * dlnR * (-z / T)
    return (n = n, dndT = dndT, dnda = dnda, h = h, dhdT = dhdT,
            I31 = n * T * h, dI31dT = dndT * T * h + n * h + n * T * dhdT, dI31da = dnda * T * h,
            aQ = mq^2 + 5 * T * h, bQ = T * h, dbQdT = h + T * dhdT, cQ = cQ, dcQdT = dcQdT)
end

_tw_is_mixture(L::HadronResonaceGas_ccbar, mq) = !hq_m2_single(L, mq)
_tw_is_mixture(::Any, mq) = false

"D_s in fm from the diffusion model (`DsT/T`, the convention of `hq_transport`)."
@inline tw_Ds(T, diffusion) = DsT(diffusion, T) / T / fmGeV

# ── options ──────────────────────────────────────────────────────────────────────────────────

"""
    TwoWayOptions

What a two-way closure carries besides `params`. Built by the closure constructors from their
keywords; `G` (a Bool type parameter) forces the generic-EOS path for a ConformalEOS.
"""
struct TwoWayOptions{G,C,BG}
    charm::C
    c::Float64            # 1 ⇔ c_M = D_s/T
    dpp_factor::Float64   # δ_ππ / τ_π
    dPP_factor::Float64   # δ_ΠΠ / τ_Π
    wv::Float64           # 2+1D corotation coupling (1 = on)
    dmn_eps::Float64      # vacuum regulator, number row only
    cM_nfloor::Float64    # optional c damping n/(n + floor); 0 = off
    nu_cap::Float64       # vacuum treatment: damp the current where |ν|/n > nu_cap; 0 = off
    nu_cap_rate::Float64  # … with relaxation rate up to (1 + nu_cap_rate)/τ_n
    α_max::Float64
    edge::Bool            # cold-edge decoupling D_s → f(T) D_s below T_edge (header, (v))
    T_edge::Float64       # f = 1 for T ≥ T_edge (GeV)
    edge_width::Float64   # f rises from edge_floor at T_edge − edge_width to 1 at T_edge (C² smoothstep)
    edge_floor::Float64   # f below the window
    corona_flux::Bool     # with edge: flux-form charm-number update below T_edge (two_way_problem only)
    background_fields::BG
end

function TwoWayOptions(; charm = nothing, cM::Bool = true, dpp_factor = 4/3, dPP_factor = 0.0,
                       vorticity::Bool = true, dmn_eps = 1e-8, cM_nfloor = 0.0, nu_cap = 0.0,
                       nu_cap_rate = 50.0,
                       α_max = 200.0, edge::Bool = true, T_edge = 0.156, edge_width = 0.020, edge_floor = 0.01,
                       corona_flux::Bool = true, background_fields = nothing, generic::Bool = false)
    (0 < edge_floor <= 1 && edge_width > 0) || error("TwoWayOptions: need 0 < edge_floor ≤ 1 and edge_width > 0")
    TwoWayOptions{generic, typeof(charm), typeof(background_fields)}(
        charm, cM ? 1.0 : 0.0, float(dpp_factor), float(dPP_factor), vorticity ? 1.0 : 0.0,
        float(dmn_eps), float(cM_nfloor), float(nu_cap), float(nu_cap_rate), float(α_max),
        edge, float(T_edge), float(edge_width), float(edge_floor), corona_flux, background_fields)
end

# conformal path? (compile-time)
@inline _tw_conf(::ConformalEOS, ::TwoWayOptions{false}) = Val(true)
@inline _tw_conf(::Any, ::TwoWayOptions) = Val(false)

# ── the medium inputs ────────────────────────────────────────────────────────────────────────

# generic: pres, dtp, dtdtp, zeta, visc, tauS, tauB, dpp, dPP
@inline function _tw_medium(T, params, o)
    th = thermodynamic(T, params.eos)
    p = th.pressure; dtp = th.pressure_derivative[1]; dtdtp = th.pressure_hessian[1]
    visc = viscosity(T, dtp, params.shear)
    tauS = τ_shear(T, dtp, params.shear)
    zeta = bulk_viscosity(T, dtp, params.bulk)
    tauB = τ_bulk(T, dtp, dtdtp, params.bulk)
    return (p, dtp, dtdtp, zeta, visc, tauS, tauB, o.dpp_factor * tauS, o.dPP_factor * tauB)
end
# conformal: aSB, visc, tauS, dpp
@inline function _tw_medium_conf(T, params, o)
    a = conformal_aSB(params.eos)
    visc = viscosity(T, 4a * T^3, params.shear)
    tauS = τ_shear(T, 4a * T^3, params.shear)
    return (a, visc, tauS, o.dpp_factor * tauS)
end

# charm inputs in the generated functions' order: n, dndT, dnda, h, dhdT, I31, dI31dT, dI31da,
# aQ, bQ, dbQdT, cQ, dcQdT, Ds   (+ the effective c)
@inline function _tw_charm(T, α, params, o)
    αs = clamp(α, -o.α_max, o.α_max)
    q = two_way_charm_thermo(T, αs, o.charm, params.diffusion.mass)
    Ds = tw_Ds(T, params.diffusion) * tw_edge_factor(T, o)
    c = o.cM_nfloor > 0 ? o.c * q.n / (q.n + o.cM_nfloor) : o.c
    return (q.n, q.dndT, q.dnda, q.h, q.dhdT, q.I31, q.dI31dT, q.dI31da,
            q.aQ, q.bQ, q.dbQdT, q.cQ, q.dcQdT, Ds), c, q
end

"""
    tw_edge_factor(T, o::TwoWayOptions) -> f

The cold-edge decoupling factor on D_s (header, item (v)): 1 for T ≥ `o.T_edge` (and whenever
`o.edge` is false), `o.edge_floor` for T ≤ T_edge − edge_width, a C² smoothstep in between.
"""
@inline function tw_edge_factor(T, o)
    (o.edge && T < o.T_edge) || return one(T)
    x = (T - (o.T_edge - o.edge_width)) / o.edge_width
    x <= 0 && return oftype(one(T), o.edge_floor)
    return o.edge_floor + (1 - o.edge_floor) * (x^3 * (10 - 15x + 6x^2))
end

# the vacuum regulator: + eps (u^τ, u^x, u^y) on the number row's own α column
@inline function _tw_reg(A::SMatrix{N,N,T}, i::Int, v) where {N,T}
    return A + SMatrix{N,N,T}(ntuple(k -> k == i + (i - 1) * N ? T(v) : zero(T), Val(N * N)))
end

# the vacuum treatment of the current: + f ν in the current row(s), with
#     w = |ν|/n  (|ν|² = ν_μν^μ, the drift speed of the charm relative to u),
#     x = max(0, w/nu_cap - 1),   f = nu_cap_rate · x²/(1 + x²).
# f ≡ 0 wherever the charm current is subluminal by the margin nu_cap (N^μ timelike), i.e. wherever the
# moment expansion can be valid; it switches on smoothly (C¹) only where |ν| ≳ n, which happens where
# n → 0 (the dilute corona / vacuum tail, Knudsen number ≳ 1) and nowhere else.
@inline function _tw_fnu(n, ν2, o)
    o.nu_cap > 0 || return zero(n)
    x = max(zero(n), sqrt(max(ν2, zero(ν2))) / (o.nu_cap * max(n, floatmin(n))) - 1)
    return o.nu_cap_rate * x^2 / (1 + x^2)
end
@inline function _tw_add(src::SVector{N,T}, i::Int, v) where {N,T}
    return src + SVector{N,T}(ntuple(k -> k == i ? T(v) : zero(T), Val(N)))
end

# ══ 0+1D Bjorken ═════════════════════════════════════════════════════════════════════════════
# Fields: S1 (T); S2 (T, π^η_η, Π) / conformal (T, π^η_η); S3 (α); S4 (α, π_Q^η_η, Π_Q);
# S5 (T, π^η_η, Π, α, π_Q^η_η, Π_Q) / conformal (T, π^η_η, α, π_Q^η_η, Π_Q).
# `background_fields(tau)` for S3/S4 returns (T, dtT).

"""
    tw_system_0p1d(sys, U, tau, params, o) -> (A_t, source)

Raw 0+1D system `sys ∈ (:ideal, :viscous, :hq_fm_bg, :hq_m2_bg, :two_way)` at state `U`.
"""
tw_system_0p1d(sys::Symbol, U, tau, params, o) = tw_system_0p1d(Val(sys), U, tau, params, o)

tw_system_0p1d(::Val{:ideal}, U, tau, params, o) = _tw_ideal0(_tw_conf(params.eos, o), U, tau, params, o)
_tw_ideal0(::Val{true}, U, tau, params, o) = ideal_0p1d_conformal_generated(U, tau, conformal_aSB(params.eos))
function _tw_ideal0(::Val{false}, U, tau, params, o)
    th = thermodynamic(U[1], params.eos)
    ideal_0p1d_generated(U, tau, th.pressure, th.pressure_derivative[1], th.pressure_hessian[1])
end

tw_system_0p1d(::Val{:viscous}, U, tau, params, o) = _tw_visc0(_tw_conf(params.eos, o), U, tau, params, o)
_tw_visc0(::Val{true}, U, tau, params, o) = viscous_0p1d_conformal_generated(U, tau, _tw_medium_conf(U[1], params, o)...)
_tw_visc0(::Val{false}, U, tau, params, o) = viscous_0p1d_generated(U, tau, _tw_medium(U[1], params, o)...)

function tw_system_0p1d(::Val{:hq_fm_bg}, U, tau, params, o)
    T, dtT = o.background_fields(tau)
    ch, _, _ = _tw_charm(T, U[1], params, o)
    At, src = hq_fm_bg_0p1d_generated(U, tau, T, dtT, ch[1:5]..., ch[14])
    return _tw_reg(At, 1, o.dmn_eps), src
end
function tw_system_0p1d(::Val{:hq_m2_bg}, U, tau, params, o)
    T, dtT = o.background_fields(tau)
    ch, _, _ = _tw_charm(T, U[1], params, o)
    At, src = hq_m2_bg_0p1d_generated(U, tau, T, dtT, ch...)
    return _tw_reg(At, 1, o.dmn_eps), src
end

tw_system_0p1d(::Val{:two_way}, U, tau, params, o) = _tw_tw0(_tw_conf(params.eos, o), U, tau, params, o)
function _tw_tw0(::Val{true}, U, tau, params, o)
    ch, c, _ = _tw_charm(U[1], U[3], params, o)
    At, src = two_way_0p1d_conformal_generated(U, tau, _tw_medium_conf(U[1], params, o)..., ch..., c)
    return _tw_reg(At, 3, o.dmn_eps), src
end
function _tw_tw0(::Val{false}, U, tau, params, o)
    ch, c, _ = _tw_charm(U[1], U[4], params, o)
    At, src = two_way_0p1d_generated(U, tau, _tw_medium(U[1], params, o)..., ch..., c)
    return _tw_reg(At, 4, o.dmn_eps), src
end

"""
    tw_bjorken_rhs(sys, params, o) -> f(U, p, tau)

Out-of-place ODE right-hand side `dU/dτ = -A_t⁻¹ source` of the 0+1D system `sys` (an
`SVector` state) — for `ODEProblem(f, U0, tspan)`.
"""
function tw_bjorken_rhs(sys::Symbol, params, o::TwoWayOptions)
    v = Val(sys)
    return (U, p, tau) -> begin
        At, src = tw_system_0p1d(v, U, tau, params, o)
        -(At \ src)
    end
end

"""
    two_way_bjorken(sys, U0, tspan; params, alg = Tsit5(), reltol = 1e-10, abstol = 1e-12,
                    charm, cM, dpp_factor, vorticity, dmn_eps, cM_nfloor, background_fields,
                    generic, solve_kwargs...)

Solve the 0+1D Bjorken system `sys ∈ (:ideal, :viscous, :hq_fm_bg, :hq_m2_bg, :two_way)` from
`U0` (field order in the comment above `tw_system_0p1d`). Returns the ODESolution.
"""
function two_way_bjorken(sys::Symbol, U0, tspan; params, alg = Tsit5(), reltol = 1e-10,
                         abstol = 1e-12, charm = nothing, cM::Bool = true, dpp_factor = 4/3,
                         dPP_factor = 0.0, dmn_eps = 1e-8, cM_nfloor = 0.0, nu_cap = 0.0,
                         nu_cap_rate = 50.0, edge::Bool = true, T_edge = 0.156, edge_width = 0.020,
                         edge_floor = 0.01, background_fields = nothing, generic::Bool = false, kwargs...)
    o = TwoWayOptions(; charm, cM, dpp_factor, dPP_factor, dmn_eps, cM_nfloor, nu_cap, nu_cap_rate,
                      edge, T_edge, edge_width, edge_floor, background_fields, generic)
    f = tw_bjorken_rhs(sys, params, o)
    prob = ODEProblem(f, SVector{length(U0)}(float.(U0)), float.(tspan))
    return solve(prob, alg; reltol, abstol, kwargs...)
end

# ══ 1+1D cylindrical ═════════════════════════════════════════════════════════════════════════
# Fields: S1 (T, u^r) = `tw_ideal_1d`; S2 `viscous_1d` (T, u^r, π^φ_φ, π^η_η, Π) / conformal
# `tw_viscous_conformal_1d` (no Π); S3 `HQ_const_BG_init` (α, ν^r); S4
# `HQ_const_BG_second_moment_init` (α, ν^r, p_l, p_φ, Π_Q); S5 `two_way_1d` (10) /
# `two_way_conformal_1d` (9). `background_fields(tau, r)` for S3/S4 returns
# (T, u^r, ∂_τT, ∂_rT, ∂_ru^r, ∂_τu^r) — the convention of the existing 1-D BG drivers.

"""
    tw_system_1p1d(sys, ϕ, tau, r, params, o) -> (A_t, A_r, source)

Raw 1+1D cylindrical system (not divided by A_t), regulator included.
"""
tw_system_1p1d(sys::Symbol, ϕ, tau, r, params, o) = tw_system_1p1d(Val(sys), ϕ, tau, r, params, o)

tw_system_1p1d(::Val{:ideal}, ϕ, tau, r, params, o) = _tw_ideal1(_tw_conf(params.eos, o), ϕ, tau, r, params, o)
_tw_ideal1(::Val{true}, ϕ, tau, r, params, o) = ideal_1p1d_conformal_generated(ϕ, tau, r, conformal_aSB(params.eos))
function _tw_ideal1(::Val{false}, ϕ, tau, r, params, o)
    th = thermodynamic(ϕ[1], params.eos)
    ideal_1p1d_generated(ϕ, tau, r, th.pressure, th.pressure_derivative[1], th.pressure_hessian[1])
end

tw_system_1p1d(::Val{:viscous}, ϕ, tau, r, params, o) = _tw_visc1(_tw_conf(params.eos, o), ϕ, tau, r, params, o)
_tw_visc1(::Val{true}, ϕ, tau, r, params, o) = viscous_1p1d_conformal_generated(ϕ, tau, r, _tw_medium_conf(ϕ[1], params, o)...)
_tw_visc1(::Val{false}, ϕ, tau, r, params, o) = viscous_1p1d_generated(ϕ, tau, r, _tw_medium(ϕ[1], params, o)...)

function tw_system_1p1d(::Val{:hq_fm_bg}, ϕ, tau, r, params, o)
    T, ur, dtT, drT, drur, dtur = o.background_fields(tau, r)
    ch, _, q = _tw_charm(T, ϕ[1], params, o)
    At, Ax, src = hq_fm_bg_1p1d_generated(ϕ, tau, r, T, ur, dtT, drT, dtur, drur, ch[1:5]..., ch[14])
    uτ = sqrt(1 + ur^2)
    return _tw_reg(At, 1, o.dmn_eps * uτ), _tw_reg(Ax, 1, o.dmn_eps * ur), _tw_add(src, 2, _tw_fnu(q.n, ϕ[2]^2 / uτ^2, o) * ϕ[2])
end
function tw_system_1p1d(::Val{:hq_m2_bg}, ϕ, tau, r, params, o)
    T, ur, dtT, drT, drur, dtur = o.background_fields(tau, r)
    ch, c, q = _tw_charm(T, ϕ[1], params, o)
    At, Ax, src = hq_m2_bg_1p1d_generated(ϕ, tau, r, T, ur, dtT, drT, dtur, drur, ch..., c)
    uτ = sqrt(1 + ur^2)
    return _tw_reg(At, 1, o.dmn_eps * uτ), _tw_reg(Ax, 1, o.dmn_eps * ur), _tw_add(src, 2, _tw_fnu(q.n, ϕ[2]^2 / uτ^2, o) * ϕ[2])
end

tw_system_1p1d(::Val{:two_way}, ϕ, tau, r, params, o) = _tw_tw1(_tw_conf(params.eos, o), ϕ, tau, r, params, o)
function _tw_tw1(::Val{true}, ϕ, tau, r, params, o)
    ch, c, q = _tw_charm(ϕ[1], ϕ[5], params, o)
    At, Ax, src = two_way_1p1d_conformal_generated(ϕ, tau, r, _tw_medium_conf(ϕ[1], params, o)..., ch..., c)
    uτ = sqrt(1 + ϕ[2]^2)
    return _tw_reg(At, 5, o.dmn_eps * uτ), _tw_reg(Ax, 5, o.dmn_eps * ϕ[2]), _tw_add(src, 6, _tw_fnu(q.n, ϕ[6]^2 / uτ^2, o) * ϕ[6])
end
function _tw_tw1(::Val{false}, ϕ, tau, r, params, o)
    ch, c, q = _tw_charm(ϕ[1], ϕ[6], params, o)
    At, Ax, src = two_way_1p1d_generated(ϕ, tau, r, _tw_medium(ϕ[1], params, o)..., ch..., c)
    uτ = sqrt(1 + ϕ[2]^2)
    return _tw_reg(At, 6, o.dmn_eps * uτ), _tw_reg(Ax, 6, o.dmn_eps * ϕ[2]), _tw_add(src, 7, _tw_fnu(q.n, ϕ[7]^2 / uτ^2, o) * ϕ[7])
end

# the closure body: A_t⁻¹ A_r, A_t⁻¹ src into Fluidum's buffers
@inline function _tw_apply1!(A_i, Source, ϕ, tau, X, params, o, v)
    At, Ax, src = tw_system_1p1d(v, ϕ, tau, X[1], params, o)
    Ainv = inv(At)
    A_mul_B!(A_i[1], Ainv, Ax)
    jgemvavx!(Source, Ainv, src)
    return nothing
end

"""
    tw_ideal_1p1d!(; generic = false)
    tw_viscous_1p1d!(; dpp_factor = 4/3, dPP_factor = 0, generic = false)
    tw_hq_fm_bg_1p1d!(; charm, background_fields, dmn_eps = 1e-8, α_max = 200)
    tw_hq_m2_bg_1p1d!(; charm, background_fields, cM = true, dmn_eps, cM_nfloor = 0, α_max)
    two_way_1p1d!(; charm, cM = true, dpp_factor = 4/3, dPP_factor = 0, dmn_eps = 1e-8,
                    cM_nfloor = 0, α_max = 200, generic = false,
                    edge = true, T_edge = 0.156, edge_width = 0.02, edge_floor = 0.01)

Closures `(A_i, Source, ϕ, tau, X, params) -> …` for `oneshoot`/`problem` (1+1D cylindrical).
`params = FluidProperties(medium_eos, shear, bulk, diffusion)`; `charm` is the charm EOS.
Field layouts: `tw_ideal_1d`, `viscous_1d` / `tw_viscous_conformal_1d`, `HQ_const_BG_init`,
`HQ_const_BG_second_moment_init`, `two_way_1d` / `two_way_conformal_1d`.
The charm closures carry the cold-edge DECOUPLING (`edge`, header (v)); for the full treatment
(with the corona number update) build the problem with `two_way_problem`.
"""
function two_way_1p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    o.charm === nothing && error("two_way_1p1d!: pass the charm EOS as `charm = …`")
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply1!(A_i, Source, ϕ, tau, X, params, o, Val(:two_way))
end
function tw_ideal_1p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply1!(A_i, Source, ϕ, tau, X, params, o, Val(:ideal))
end
function tw_viscous_1p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply1!(A_i, Source, ϕ, tau, X, params, o, Val(:viscous))
end
function tw_hq_fm_bg_1p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (o.charm === nothing || o.background_fields === nothing) && error("tw_hq_fm_bg_1p1d!: needs `charm` and `background_fields`")
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply1!(A_i, Source, ϕ, tau, X, params, o, Val(:hq_fm_bg))
end
function tw_hq_m2_bg_1p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (o.charm === nothing || o.background_fields === nothing) && error("tw_hq_m2_bg_1p1d!: needs `charm` and `background_fields`")
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply1!(A_i, Source, ϕ, tau, X, params, o, Val(:hq_m2_bg))
end

# ══ 2+1D Cartesian ═══════════════════════════════════════════════════════════════════════════
# Fields: S1 `tw_ideal_2d` (T, u^x, u^y); S2 `viscous_2d` (T, u^x, u^y, π^{yy}, π^η_η, π^{xy}, Π) /
# `tw_viscous_conformal_2d` (no Π); S3 `HQ_2p1d_BG_init` (α, ν^x, ν^y); S4 `HQ_2p1d_BG_m2_7f_init`
# (α, ν^x, ν^y, π_Q^{xx}, π_Q^{xy}, π_Q^{yy}, Π_Q); S5 `two_way_2d` (14) / `two_way_conformal_2d` (13).
# `background_fields(tau, x, y)` for S3/S4 returns (T, u^x, u^y, ∂_τT, ∂_xT, ∂_yT, ∂_τu^x, ∂_xu^x,
# ∂_yu^x, ∂_τu^y, ∂_xu^y, ∂_yu^y) — the convention of the existing 2-D BG drivers.

"""
    tw_system_2p1d(sys, ϕ, tau, x, y, params, o) -> (A_t, A_x, A_y, source)

Raw 2+1D Cartesian system (not divided by A_t), regulator included.
"""
tw_system_2p1d(sys::Symbol, ϕ, tau, x, y, params, o) = tw_system_2p1d(Val(sys), ϕ, tau, x, y, params, o)

tw_system_2p1d(::Val{:ideal}, ϕ, tau, x, y, params, o) = _tw_ideal2(_tw_conf(params.eos, o), ϕ, tau, params, o)
_tw_ideal2(::Val{true}, ϕ, tau, params, o) = ideal_2p1d_conformal_generated(ϕ, tau, conformal_aSB(params.eos))
function _tw_ideal2(::Val{false}, ϕ, tau, params, o)
    th = thermodynamic(ϕ[1], params.eos)
    ideal_2p1d_generated(ϕ, tau, th.pressure, th.pressure_derivative[1], th.pressure_hessian[1])
end

tw_system_2p1d(::Val{:viscous}, ϕ, tau, x, y, params, o) = _tw_visc2(_tw_conf(params.eos, o), ϕ, tau, params, o)
_tw_visc2(::Val{true}, ϕ, tau, params, o) = viscous_2p1d_conformal_generated(ϕ, tau, _tw_medium_conf(ϕ[1], params, o)...)
_tw_visc2(::Val{false}, ϕ, tau, params, o) = viscous_2p1d_generated(ϕ, tau, _tw_medium(ϕ[1], params, o)...)

function tw_system_2p1d(::Val{:hq_fm_bg}, ϕ, tau, x, y, params, o)
    bg = o.background_fields(tau, x, y)
    T, ux, uy = bg[1], bg[2], bg[3]
    ch, _, q = _tw_charm(T, ϕ[1], params, o)
    At, Ax, Ay, src = hq_fm_bg_2p1d_generated(ϕ, tau, bg..., ch[1:5]..., ch[14])
    uτ = sqrt(1 + ux^2 + uy^2); e = o.dmn_eps
    f = _tw_fnu(q.n, ϕ[2]^2 + ϕ[3]^2 - ((ux * ϕ[2] + uy * ϕ[3]) / uτ)^2, o)
    return _tw_reg(At, 1, e * uτ), _tw_reg(Ax, 1, e * ux), _tw_reg(Ay, 1, e * uy),
           _tw_add(_tw_add(src, 2, f * ϕ[2]), 3, f * ϕ[3])
end
function tw_system_2p1d(::Val{:hq_m2_bg}, ϕ, tau, x, y, params, o)
    bg = o.background_fields(tau, x, y)
    T, ux, uy = bg[1], bg[2], bg[3]
    ch, c, q = _tw_charm(T, ϕ[1], params, o)
    At, Ax, Ay, src = hq_m2_bg_2p1d_generated(ϕ, tau, bg..., ch..., c, o.wv)
    uτ = sqrt(1 + ux^2 + uy^2); e = o.dmn_eps
    f = _tw_fnu(q.n, ϕ[2]^2 + ϕ[3]^2 - ((ux * ϕ[2] + uy * ϕ[3]) / uτ)^2, o)
    return _tw_reg(At, 1, e * uτ), _tw_reg(Ax, 1, e * ux), _tw_reg(Ay, 1, e * uy),
           _tw_add(_tw_add(src, 2, f * ϕ[2]), 3, f * ϕ[3])
end

tw_system_2p1d(::Val{:two_way}, ϕ, tau, x, y, params, o) = _tw_tw2(_tw_conf(params.eos, o), ϕ, tau, params, o)
function _tw_tw2(::Val{true}, ϕ, tau, params, o)
    ch, c, q = _tw_charm(ϕ[1], ϕ[7], params, o)
    At, Ax, Ay, src = two_way_2p1d_conformal_generated(ϕ, tau, _tw_medium_conf(ϕ[1], params, o)..., ch..., c, o.wv)
    uτ = sqrt(1 + ϕ[2]^2 + ϕ[3]^2); e = o.dmn_eps
    f = _tw_fnu(q.n, ϕ[8]^2 + ϕ[9]^2 - ((ϕ[2] * ϕ[8] + ϕ[3] * ϕ[9]) / uτ)^2, o)
    return _tw_reg(At, 7, e * uτ), _tw_reg(Ax, 7, e * ϕ[2]), _tw_reg(Ay, 7, e * ϕ[3]),
           _tw_add(_tw_add(src, 8, f * ϕ[8]), 9, f * ϕ[9])
end
function _tw_tw2(::Val{false}, ϕ, tau, params, o)
    ch, c, q = _tw_charm(ϕ[1], ϕ[8], params, o)
    At, Ax, Ay, src = two_way_2p1d_generated(ϕ, tau, _tw_medium(ϕ[1], params, o)..., ch..., c, o.wv)
    uτ = sqrt(1 + ϕ[2]^2 + ϕ[3]^2); e = o.dmn_eps
    f = _tw_fnu(q.n, ϕ[9]^2 + ϕ[10]^2 - ((ϕ[2] * ϕ[9] + ϕ[3] * ϕ[10]) / uτ)^2, o)
    return _tw_reg(At, 8, e * uτ), _tw_reg(Ax, 8, e * ϕ[2]), _tw_reg(Ay, 8, e * ϕ[3]),
           _tw_add(_tw_add(src, 9, f * ϕ[9]), 10, f * ϕ[10])
end

@inline function _tw_apply2!(A_i, Source, ϕ, tau, X, params, o, v)
    At, Ax, Ay, src = tw_system_2p1d(v, ϕ, tau, X[1], X[2], params, o)
    Ainv = inv(At)
    A_mul_B!(A_i[1], Ainv, Ax)
    A_mul_B!(A_i[2], Ainv, Ay)
    jgemvavx!(Source, Ainv, src)
    return nothing
end

"""
    tw_ideal_2p1d!(; generic = false)
    tw_viscous_2p1d!(; dpp_factor = 4/3, dPP_factor = 0, generic = false)
    tw_hq_fm_bg_2p1d!(; charm, background_fields, dmn_eps = 1e-8, α_max = 200)
    tw_hq_m2_bg_2p1d!(; charm, background_fields, cM = true, vorticity = true, …)
    two_way_2p1d!(; charm, cM = true, dpp_factor = 4/3, vorticity = true, dmn_eps = 1e-8,
                    cM_nfloor = 0, α_max = 200, generic = false)

The 2+1D Cartesian closures (see the 1+1D docstring; `vorticity` is the corotation coupling
of the charm second moment, on by default). Layouts: `tw_ideal_2d`,
`viscous_2d` / `tw_viscous_conformal_2d`, `HQ_2p1d_BG_init`, `HQ_2p1d_BG_m2_7f_init`,
`two_way_2d` / `two_way_conformal_2d`.
"""
function two_way_2p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    o.charm === nothing && error("two_way_2p1d!: pass the charm EOS as `charm = …`")
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply2!(A_i, Source, ϕ, tau, X, params, o, Val(:two_way))
end
function tw_ideal_2p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply2!(A_i, Source, ϕ, tau, X, params, o, Val(:ideal))
end
function tw_viscous_2p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply2!(A_i, Source, ϕ, tau, X, params, o, Val(:viscous))
end
function tw_hq_fm_bg_2p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (o.charm === nothing || o.background_fields === nothing) && error("tw_hq_fm_bg_2p1d!: needs `charm` and `background_fields`")
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply2!(A_i, Source, ϕ, tau, X, params, o, Val(:hq_fm_bg))
end
function tw_hq_m2_bg_2p1d!(; kw...)
    o = TwoWayOptions(; kw...)
    (o.charm === nothing || o.background_fields === nothing) && error("tw_hq_m2_bg_2p1d!: needs `charm` and `background_fields`")
    (A_i, Source, ϕ, tau, X, params) -> _tw_apply2!(A_i, Source, ϕ, tau, X, params, o, Val(:hq_m2_bg))
end

# ══ the corona number update and `two_way_problem` ═══════════════════════════════════════════
# Below T_edge the quasi-linear α update of discreteFields.jl is replaced (blended over the window
# with weight w = 1 − S((T − T_edge + edge_width)/edge_width), S the smoothstep of tw_edge_factor) by
# the FLUX-FORM update of the charm number,
#     d(τ r N^τ)_i/dτ = −(G_{i+1/2} − G_{i−1/2})/Δr,   G = τ r [ (N^r_L + N^r_R)/2 − a (N^τ_R − N^τ_L)/2 ],
#     a = min(1, max |N^r|/N^τ over the two cells) (local Lax–Friedrichs: positivity-preserving for the number),
# solved for α̇ given the scheme's own Ṫ, u̇, ν̇ (∂_αn regulated by dmn_eps, as in the number row).
# (2+1D Cartesian: the same with τ N^x, τ N^y on the four faces.) The quasi-linear update is not
# positivity-preserving at an under-resolved charm edge. A cell next to a neighbour with a much larger
# current is charged half of that neighbour's ∂ν (central difference), its n collapses (α → −10³), and
# the scheme's dissipation, which acts on α, then drains the next cell. For T ≥ T_edge nothing changes.
# The fluid rows keep the quasi-linear α̇, so the charm energy ∂ε_Q/∂α in the medium energy row is off
# by O(Δα̇) in the corona only.

"n, ∂_Tn, ∂_αn of the charm EOS (no h): the closed form or the list's `thermodynamic`"
@inline function _tw_ndens(T, α, x::BoltzmannCharm)
    n, dn, _, _ = _tw_species(T, α, x.m, x.g, 1)
    return n, dn, n
end
@inline function _tw_ndens(T, α, x::HadronResonaceGas_ccbar)
    th = thermodynamic(T, α, x)
    return th.pressure, th.pressure_derivative[1], th.pressure_derivative[2]
end

"the corona weight: 0 for T ≥ T_edge, 1 for T ≤ T_edge − edge_width, smoothstep in between"
@inline function _tw_corona_weight(T, o)
    (o.edge && o.corona_flux && T < o.T_edge) || return zero(T)
    x = (T - (o.T_edge - o.edge_width)) / o.edge_width
    x <= 0 && return one(T)
    return 1 - x^3 * (10 - 15x + 6x^2)
end

# a = the larger drift |N^r|/N^τ of the two cells, capped at 1 (positivity holds wherever the charm
# drift is subluminal; beyond that the moment description is void, and an uncapped a — |ν| ≫ n in a
# near-vacuum cell — makes the dissipation stiff enough to blow up an explicit step)
@inline function _tw_llf(UL, UR, FL, FR)
    a = min(one(UL), max(UL > 0 ? abs(FL) / UL : one(UL), UR > 0 ? abs(FR) / UR : one(UR)))
    return (FL + FR) / 2 - a * (UR - UL) / 2
end

# N^τ, N^r of 1+1D cell i
@inline function _tw_N1(u, i, l, o, iα, iν)
    ur = u[2, i]; uτ = sqrt(1 + ur^2)
    n, _, _ = _tw_ndens(u[1, i], clamp(u[iα, i], -o.α_max, o.α_max), l)
    return n * uτ + ur * u[iν, i] / uτ, n * ur + u[iν, i]
end
# N^τ, N^x, N^y of 2+1D cell (i, j)
@inline function _tw_N2(u, i, j, l, o, iα, iν)
    ux = u[2, i, j]; uy = u[3, i, j]; uτ = sqrt(1 + ux^2 + uy^2); νx = u[iν, i, j]; νy = u[iν+1, i, j]
    n, _, _ = _tw_ndens(u[1, i, j], clamp(u[iα, i, j], -o.α_max, o.α_max), l)
    return n * uτ + (ux * νx + uy * νy) / uτ, n * ux + νx, n * uy + νy
end

# α̇ that gives the target d(τ·vol·N^τ)/dτ (vol = r in 1+1D, 1 in 2+1D), cell state and rates given
@inline function _tw_alpha_rate(target, vol, τ, T, α, uτ, u̇τ, dντ, Ṫ, Nt, l, o)
    n, dndT, dnda = _tw_ndens(T, clamp(α, -o.α_max, o.α_max), l)
    dna = abs(α) < o.α_max ? dnda : zero(dnda)
    return ((target - vol * Nt) / (τ * vol) - dndT * Ṫ * uτ - n * u̇τ - dντ) / (uτ * (dna + o.dmn_eps))
end

"""
    TwoWayCoronaRHS

The right-hand side of `two_way_problem`: Fluidum's quasi-linear upwinding (`f`) followed by the
corona number update (comment above). Callable `(du, u, p, t)`; allocation-free.
"""
struct TwoWayCoronaRHS{FN,D,O,R}
    f::FN
    d::D
    o::O
    rr::R                 # 1+1D: cell-centre radii (ghosted); 2+1D: nothing
    dx::Float64
    dy::Float64
    iα::Int
    iν::Int               # ν^r (1+1D) or ν^x (2+1D; ν^y = iν + 1)
    sd::Int
end

function (c::TwoWayCoronaRHS)(du, u, p, t)
    c.f(du, u, p, t)
    c.sd == 1 ? _tw_corona_1d!(du, u, t, c) : _tw_corona_2d!(du, u, t, c)
    boundary_condition!(du, c.d)
    return nothing
end

function _tw_corona_1d!(du, u, τ, c)
    o = c.o; l = _tw_charm_list(o.charm); rr = c.rr; dr = c.dx; iα = c.iα; iν = c.iν
    N = size(u, 2)
    @inbounds for i in 2:N-1
        w = _tw_corona_weight(u[1, i], o)
        iszero(w) && continue
        NtC, NrC = _tw_N1(u, i, l, o, iα, iν)
        rm = rr[i] - dr / 2; rp = rr[i] + dr / 2
        Gm = if rm <= 0
            zero(NtC)
        else
            NtL, NrL = _tw_N1(u, i - 1, l, o, iα, iν); τ * rm * _tw_llf(NtL, NtC, NrL, NrC)
        end
        NtR, NrR = _tw_N1(u, i + 1, l, o, iα, iν)
        Gp = τ * rp * _tw_llf(NtC, NtR, NrC, NrR)
        ur = u[2, i]; uτ = sqrt(1 + ur^2); ν = u[iν, i]
        u̇ = du[2, i]; ν̇ = du[iν, i]; u̇τ = ur * u̇ / uτ
        dντ = (u̇ * ν + ur * ν̇) / uτ - ur * ν * u̇τ / uτ^2
        αff = _tw_alpha_rate(-(Gp - Gm) / dr, rr[i], τ, u[1, i], u[iα, i], uτ, u̇τ, dντ, du[1, i], NtC, l, o)
        du[iα, i] = (1 - w) * du[iα, i] + w * αff
    end
    return nothing
end

function _tw_corona_2d!(du, u, τ, c)
    o = c.o; l = _tw_charm_list(o.charm); iα = c.iα; iν = c.iν
    Nx, Ny = size(u, 2), size(u, 3)
    @inbounds for j in 2:Ny-1, i in 2:Nx-1
        w = _tw_corona_weight(u[1, i, j], o)
        iszero(w) && continue
        UC, XC, YC = _tw_N2(u, i, j, l, o, iα, iν)
        UL, XL, _ = _tw_N2(u, i - 1, j, l, o, iα, iν); UR, XR, _ = _tw_N2(u, i + 1, j, l, o, iα, iν)
        UD, _, YD = _tw_N2(u, i, j - 1, l, o, iα, iν); UU, _, YU = _tw_N2(u, i, j + 1, l, o, iα, iν)
        target = -τ * ((_tw_llf(UC, UR, XC, XR) - _tw_llf(UL, UC, XL, XC)) / c.dx +
                       (_tw_llf(UC, UU, YC, YU) - _tw_llf(UD, UC, YD, YC)) / c.dy)
        ux = u[2, i, j]; uy = u[3, i, j]; uτ = sqrt(1 + ux^2 + uy^2); νx = u[iν, i, j]; νy = u[iν+1, i, j]
        u̇x = du[2, i, j]; u̇y = du[3, i, j]
        u̇τ = (ux * u̇x + uy * u̇y) / uτ; ντ = (ux * νx + uy * νy) / uτ
        dντ = (u̇x * νx + ux * du[iν, i, j] + u̇y * νy + uy * du[iν+1, i, j]) / uτ - ντ * u̇τ / uτ
        αff = _tw_alpha_rate(target, one(τ), τ, u[1, i, j], u[iα, i, j], uτ, u̇τ, dντ, du[1, i, j], UC, l, o)
        du[iα, i, j] = (1 - w) * du[iα, i, j] + w * αff
    end
    return nothing
end

"""
    two_way_problem(d::DiscreteFields, params, init, tspan; charm, cM = true, kw...) -> ODEProblem

The two-way system (S5) on `d` (layout `two_way_conformal_1d` / `two_way_1d` in 1+1D cylindrical,
`two_way_conformal_2d` / `two_way_2d` in 2+1D Cartesian) WITH the full cold-edge treatment: the
decoupling D_s → f(T) D_s below T_edge (in the closure, as for `two_way_1p1d!` / `two_way_2p1d!`)
AND the flux-form corona number update (`TwoWayCoronaRHS`). `kw` are the `TwoWayOptions` keywords;
`edge = false` gives the untreated system (= `problem(d, two_way_1p1d!(; edge = false, …), …)`),
`corona_flux = false` the decoupling alone. Solve with any OrdinaryDiffEq algorithm
(e.g. Tsit5, reltol 1e-8, abstol 1e-10, dtmax 0.01).
"""
function two_way_problem(d::DiscreteFields, params, init, tspan; kw...)
    o = TwoWayOptions(; kw...)
    o.charm === nothing && error("two_way_problem: pass the charm EOS as `charm = …`")
    X = d.discretization
    dxs = Δx(X)
    sd = length(dxs)
    conf = _tw_conf(params.eos, o) isa Val{true}
    m = sd == 1 ? ((A_i, Source, ϕ, tau, Xc, prm) -> _tw_apply1!(A_i, Source, ϕ, tau, Xc, prm, o, Val(:two_way))) :
                  ((A_i, Source, ϕ, tau, Xc, prm) -> _tw_apply2!(A_i, Source, ϕ, tau, Xc, prm, o, Val(:two_way)))
    base = problem(d, m, params, init, tspan)
    (o.edge && o.corona_flux) || return base
    nf = size(init, 1)
    if sd == 1
        nf == (conf ? 9 : 10) || error("two_way_problem: 1+1D layout has $nf fields, expected $(conf ? 9 : 10)")
        iα = conf ? 5 : 6
        rr = [X.grid[i][1] for i in eachindex(X.grid)]
        g = TwoWayCoronaRHS(base.f, d, o, rr, float(dxs[1]), float(dxs[1]), iα, iα + 1, 1)
    elseif sd == 2
        nf == (conf ? 13 : 14) || error("two_way_problem: 2+1D layout has $nf fields, expected $(conf ? 13 : 14)")
        iα = conf ? 7 : 8
        g = TwoWayCoronaRHS(base.f, d, o, nothing, float(dxs[1]), float(dxs[2]), iα, iα + 1, 2)
    else
        error("two_way_problem: 1+1D or 2+1D only")
    end
    return ODEProblem(g, base.u0, base.tspan)
end

# ══ initial conditions ═══════════════════════════════════════════════════════════════════════

"""
    tw_alpha_from_density(T, n, charm; α_min = -150.0) -> α

The fugacity at which the charm EOS has density `n` at temperature `T` (unit-charge lists:
α = log(n / n(T, 0)); with a doubly charmed state, a few Newton steps). Clamped
below at `α_min` so a vacuum tail stays finite.
"""
function tw_alpha_from_density(T, n, charm; α_min = -150.0)
    n <= 0 && return α_min
    l = _tw_charm_list(charm)
    α = log(n / _tw_charm_core(T, 0.0, l).n)
    for _ in 1:30
        c = _tw_charm_core(T, α, l)
        dα = (c.n - n) / c.dnda
        α -= dα
        abs(dα) < 1e-14 * max(1, abs(α)) && break
    end
    return max(α, α_min)
end
