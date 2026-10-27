# This branch on top of the upstream code: what phd-git uses and upstream does not carry. Every
# function here is a new name or a new method; nothing defined elsewhere in the package is changed.
# Included last from Fluidum.jl.

# ── Process-wide switches ────────────────────────────────────────────────────────────────────────
# Each is a Ref set in __init__ from the environment, so it is read once per process at load time.
# Set the env var before `using Fluidum`, or assign the Ref before building the matrix.
#
# HQ_TAUN_SCALE: multiplier on the heavy-quark diffusion relaxation time τ_n, applied by `hq_taun` /
# `hq_transport`. Default 1.0 is the bare τ_n of Eq. (30) of 2205.07692. Diagnostic only.
const HQ_TAUN_SCALE = Ref(1.0)
# FLUIDUM_HQ_MOMENTUM_DIMS: 3 (default) or 2. With 2, charm τ_n/τ_M are matched in the 2-D momentum
# measure of a transverse-only Langevin ensemble (see `hq_transport`). Diagnostic.
const HQ_MOMENTUM_DIMS = Ref(3)
# FLUIDUM_HQ_M2_MIXTURE: `hq_m2_mixture` uses the species sums of a multi-species charm list (true,
# the default, env "1") or the single-quark forms (env "0"). Diagnostic.
const HQ_M2_MIXTURE = Ref(true)

function _env_bool(name, default::String)
    s = get(ENV, name, default)
    s == "0" && return false
    s == "1" && return true
    error("$name must be \"0\" or \"1\", got \"$s\"")
end

function __init__()
    HQ_TAUN_SCALE[] = parse(Float64, get(ENV, "HQ_TAUN_SCALE", "1.0"))
    HQ_MOMENTUM_DIMS[] = parse(Int, get(ENV, "FLUIDUM_HQ_MOMENTUM_DIMS", "3"))
    HQ_MOMENTUM_DIMS[] in (2, 3) || error("FLUIDUM_HQ_MOMENTUM_DIMS must be 2 or 3")
    HQ_MOMENTUM_DIMS[] == 3 || @warn "Fluidum: charm τ_n/τ_M matched in the 2-D momentum measure, diagnostic mode, not production"
    HQ_TAUN_SCALE[] == 1.0 || @warn "Fluidum: τ_n scaled by HQ_TAUN_SCALE, diagnostic mode, not production" scale=HQ_TAUN_SCALE[]
    HQ_M2_MIXTURE[] = _env_bool("FLUIDUM_HQ_M2_MIXTURE", "1")
    HQ_M2_MIXTURE[] || @warn "Fluidum: charm second moment with single-quark coefficients for multi-species lists (FLUIDUM_HQ_M2_MIXTURE=0) (diagnostic)"
end

# ── Equation of state and transport ──────────────────────────────────────────────────────────────
@inline @fastmath pressure(T, x::Heavy_Quark) = pressure(T, x.eos)
@inline pressure_derivative(T, v::Val{1}, x::Heavy_Quark) = pressure_derivative(T, v, x.eos)
@inline pressure_derivative(T, v::Val{2}, x::Heavy_Quark) = pressure_derivative(T, v, x.eos)

function thermodynamic_perturbation(T::N, eos::ConformalEOS{S}) where {N,S}
    TT = float(promote_type(N,S))
    a = TT(a_SB(eos)) * TT(fmGeV^3)
    ThermodynamicPerturbation(a*T^4, (4a*T^3,), (12a*T^2,), (24a*T,))
end

# α-independent: (P, (∂_T P, ∂_α P), (∂²_T P, ∂_T∂_α P, ∂²_α P)) with every α-derivative zero.
function thermodynamic(T::N, α::M, eos::ConformalEOS{S}) where {N,M,S}
    th = thermodynamic(T, eos)
    TT = float(promote_type(N,M,S))
    z = zero(TT)
    Thermodynamic{TT,2,3}(TT(th.pressure), (TT(th.pressure_derivative[1]), z),
                          (TT(th.pressure_hessian[1]), z, z))
end

# The earlier constructor name. The linear law is frozen below T_fo = 0.156 GeV, the LP1 transport
# contract D_s T = 1.765 max(T, T_fo) - 0.159; upstream's LinearDiffusion freezes at 0.180 by default.
HQdiffusion(DsT::T, mass::T) where {T<:Real} = ConstDiffusion(DsT, mass)
HQdiffusion(slope::T, offset::T, mass::T) where {T<:Real} = LinearDiffusion(slope, offset, mass; Tc = 0.156)

# ── Charm relaxation times ───────────────────────────────────────────────────────────────────────
# phd-git's drivers take τ_n, τ_M, η_M, h from `hq_taun` / `hq_transport`, which apply HQ_TAUN_SCALE
# and FLUIDUM_HQ_MOMENTUM_DIMS. The package's charm matrices use `charm_moments` and see neither.
#
# FLUIDUM_HQ_MOMENTUM_DIMS=2 matches the charm sector in the 2-D momentum measure, for comparison with
# a Langevin ensemble that evolves transverse momenta only. Its equilibrium is the 2-D Jüttner
# p dp e^{-E/T}, and its ℓ-th moment relaxes at λ_ℓ(2D)·η_D instead of the 3-D λ_ℓ = K_{ℓ+1}/K_{ℓ+2}
# encoded in τ_n = D_s z K₃/K₂ and τ_M = D_s z K₄/(2K₃). In that measure every moment
# ∫_M^∞ E^n e^{-E/T} dE is an incomplete Gamma of integer order, so
#     λ₁(2D) = z(z+1)/(z²+3z+3),      λ₂(2D) = z(z²+3z+3)/(z³+6z²+15z+15),
# and τ_n^{2D} = τ_n·(K₂/K₃)/λ₁(2D),  τ_M^{2D} = τ_M·(K₃/K₄)/λ₂(2D).
# η_M stays tied to τ_n (η_M = T τ_n/2); κ = D_s n is dimension-independent. The tensor structure of
# the matrices stays three-dimensional. Diagnostic only.

"""
    hq_taun(T, α, eos, diffusion)

The charm diffusion relaxation time τ_n (fm) every charm matrix uses: the bare Eq. (30) τ_n of
2205.07692 (`τ_diffusion_hadron`) times `HQ_TAUN_SCALE[]`, rescaled to the 2-D momentum measure
when `HQ_MOMENTUM_DIMS[] == 2`. Equal to `hq_transport(...).taun`.
"""
function hq_taun(T, α, eos, diffusion)
    taun = τ_diffusion_hadron(T, α, eos, diffusion) * HQ_TAUN_SCALE[]
    if HQ_MOMENTUM_DIMS[] == 2
        z = diffusion.mass / T
        taun *= (besselk(2, z) / besselk(3, z)) / (z * (z + 1) / (z^2 + 3z + 3))
    end
    return taun
end

"""
    hq_transport(T, α, eos, diffusion; use_NR_tauM = false) -> (; taun, tauM, etaM, h, Ds)

The charm transport coefficients of the first and second moment, tied to one τ_n:

    τ_M = τ_n K₄K₂/(2K₃²)   (use_NR_tauM: τ_n/2),    η_M = T τ_n/2,
    D_s = DsT/T (fm),        h = τ_n T/D_s           (so τ_n = D_s h/T exactly),

with `HQ_TAUN_SCALE` and `FLUIDUM_HQ_MOMENTUM_DIMS` applied (z = m/T, m = `diffusion.mass`). The
Fokker–Planck to hydro matching fixes the ratio τ_M/τ_n = ½K₄K₂/K₃² (≈ 0.550 at z = 9.615); tying
τ_M to τ_n keeps that ratio for any normalisation of τ_n.
"""
function hq_transport(T, α, eos, diffusion; use_NR_tauM::Bool = false)
    taun = τ_diffusion_hadron(T, α, eos, diffusion) * HQ_TAUN_SCALE[]
    z  = diffusion.mass / T
    K2 = besselk(2, z); K3 = besselk(3, z); K4 = besselk(4, z)
    tauM = use_NR_tauM ? taun / 2 : taun * K4 * K2 / (2 * K3^2)
    if HQ_MOMENTUM_DIMS[] == 2
        taun *= (K2 / K3) / (z * (z + 1) / (z^2 + 3z + 3))
        tauM *= (K3 / K4) / (z * (z^2 + 3z + 3) / (z^3 + 6z^2 + 15z + 15))
    end
    etaM = T * taun / 2
    Ds = DsT(diffusion, T) / T / fmGeV
    return (taun = taun, tauM = tauM, etaM = etaM, h = taun * T / Ds, Ds = Ds)
end


# ── EOS coefficients of the consistent second moment for a charm list ────────────────────────────
# Besides h, the consistent second moment needs: I31 (its T-log-derivative A/B multiplies D ln T in
# the trace drive), A = I41/P0 (the a·ν completion terms) and C = I52/I42 (τ_M = D_s C/(2T) and
# D ln C). They are derived for one Boltzmann species. For a list of species i (mass m_i, degeneracy
# g_i, charm charge q_i) at a common (T, α) with one ν and one δM field they are species sums, with
# n = Σ q n_i and δf_i ∝ q_i f0_i in the ν and π channels, so h = Σ q² n_i h_i / Σ q² n_i:
#     I31 = n T h          ⇒ A/B = ∂ln I31/∂ln T|_α = T ∂_T n/n + 1 + T ∂_T h/h,   ∂_α ln I31 = (Σ q² n_i)/n
#     A   = Σ q² n_i (m_i² + 5T h_i) / Σ q² n_i  =  ⟨m²⟩ + 5T h,         ⟨m²⟩ = Σ q² n_i m_i² / Σ q² n_i
#     C   = Σ q² I42_i C_i / Σ q² I42_i,   I42_i = T² n_i h_i,   C_i = m_i K4(z_i)/K3(z_i)
# with every T-derivative analytic at fixed α. For one species of mass m these reduce to
# A/B = 5 + zK2/K3, A = m² + 5Th, C = m K4/K3.
# The derivation assumes unit charge. A q = 2 state (J/ψ) enters with the weights above; a
# charge-correct second moment would also carry ∂_α C, ∂_α h, and Σ q I31_i instead of nTh.

"""
    with_m2_mixture(f, on::Bool)

Run `f()` with `HQ_M2_MIXTURE[] = on` and restore it afterwards (closures read the Ref when they are
built, so build them inside `f`). `on = false` uses the single-quark coefficients.
"""
function with_m2_mixture(f, on::Bool)
    old = HQ_M2_MIXTURE[]
    HQ_M2_MIXTURE[] = on
    try
        return f()
    finally
        HQ_M2_MIXTURE[] = old
    end
end

"""
    hq_m2_single(list, mq) -> Bool

True when the charm `list` is one unit-charge species of mass `mq`. The second-moment coefficients
are then the single-quark closed forms.
"""
function hq_m2_single(list::HadronResonaceGas_ccbar, mq)
    length(list.particle_list) == 1 || return false
    p = list.particle_list[1]
    return round(p.Nc + p.Nac) == 1 && p.Mass == mq
end

"""
    hq_m2_mixture(T, α, list, mq) -> NamedTuple

The consistent second moment's EOS coefficients of a charm `list` (a `HadronResonaceGas_ccbar`) at
(T, α), from the list's species sums (comment block above):
`h, dhdT` (the q²n-weighted enthalpy per particle, GeV, and its ∂_T at fixed α),
`m2bar` = ⟨m²⟩ (A = m2bar + 5Th), `C`, `dCdT`, `dlnC_dlnT`, `ABr` = ∂ln I31/∂ln T|_α,
`aI31` = ∂_α ln I31 (= Σq²n_i/n; 1 for unit charges), and `Cratio` = C / (h K4K2/K3²|_{z = mq/T}),
the factor that turns `hq_transport`'s single-quark-ratio τ_M into the list's τ_M = D_s C/(2T)
(keeping the diagnostic knobs' factors). Used by the m2 drivers and the two-way wrappers. Does not
allocate.
"""
function hq_m2_mixture(T, α, list::HadronResonaceGas_ccbar, mq)
    R = float(promote_type(typeof(T), typeof(α)))
    # common scale e^{-sh}, sh = max_i (q_i α − z_i): every weight is a ratio, so nothing underflows
    # (deep vacuum α → −∞, cold T) and no species is dropped
    sh = typemin(R)
    for p in list.particle_list
        sh = max(sh, round(p.Nc + p.Nac) * α - p.Mass / T)
    end
    N1 = zero(R); dN1 = zero(R); D = zero(R); dD = zero(R)
    Hn = zero(R); dHn = zero(R); M2 = zero(R); WC = zero(R); dWC = zero(R)
    for p in list.particle_list
        q = round(p.Nc + p.Nac); m = p.Mass; g = p.Degeneracy
        z = m / T
        K1 = besselk1x(z); K2 = besselkx(2, z)          # exponentially scaled; the e^z cancels in ratios
        K3 = K1 + 4 / z * K2
        K4 = K2 + 6 / z * K3
        K5 = K3 + 8 / z * K4
        ex = exp(q * α - z - sh)
        ni  = g * m^2 * T * K2 * ex / (2π^2) * fmGeV3                        # fm⁻³ × e^{-sh}
        dni = g * ex * m^2 * (m * K1 + 2 * T * K2 + m * K3) / (4π^2 * T) * fmGeV3   # ∂_T at fixed α
        hi  = m * K3 / K2
        dhi = -(z / T) * m * ((-(K2 + K4) / 2) * K2 - K3 * (-(K1 + K3) / 2)) / K2^2
        Ci  = m * K4 / K3
        dCi = -(z / T) * m * ((-(K3 + K5) / 2) * K3 - K4 * (-(K2 + K4) / 2)) / K3^2
        w = q^2 * ni; dw = q^2 * dni
        N1 += q * ni; dN1 += q * dni
        D += w; dD += dw
        Hn += w * hi; dHn += dw * hi + w * dhi
        M2 += w * m^2
        WC += w * hi * Ci; dWC += (dw * hi + w * dhi) * Ci + w * hi * dCi
    end
    h = Hn / D; dhdT = (dHn - h * dD) / D
    C = WC / Hn; dCdT = (dWC - C * dHn) / Hn
    zq = mq / T
    k1 = besselk1x(zq); k2 = besselkx(2, zq); k3 = k1 + 4 / zq * k2; k4 = k2 + 6 / zq * k3
    return (h = h, dhdT = dhdT, m2bar = M2 / D,
            C = C, dCdT = dCdT, dlnC_dlnT = T * dCdT / C,
            ABr = T * (dN1 * T * h + N1 * h + N1 * T * dhdT) / (N1 * T * h),
            aI31 = D / N1, Cratio = C / (h * k4 * k2 / k3^2))
end


function normalization_without_degeneracy(T,μ,x::Heavy_Quark)
    norm = zero(float(promote_type(typeof(T), typeof(μ))))
    for i in x.hadron_list.particle_list
        m = i.Mass
        q = i.Nc + i.Nac
        n = free_hadron(T,μ,1,q; m = m)[1]
        norm += q^2*n
    end
    return norm;
end


# The ratio normalization / normalization_without_degeneracy (a ÷g factor). No matrix uses it.
second_moment_transport_normalization(T,μ,x) = one(promote_type(typeof(T), typeof(μ)))

function second_moment_transport_normalization(T,μ,x::Heavy_Quark)
    bare_norm = normalization_without_degeneracy(T,μ,x)
    bare_norm == 0 && return one(promote_type(typeof(T), typeof(μ)))
    return normalization(T,μ,x) / bare_norm
end

# ── Matrix closures with the keyword the other modes take ────────────────────────────────────────
# `dmn_eps` is accepted and ignored: the ideal fugacity matrix has no regulator on dn/dα, and the
# density matrix has no dn/dα at all.
function matrix1d_ideal_HQ_fugacity!(;dmn_eps=1e-6)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_ideal_HQ_fugacity!(A_i, Source, ϕ, tau, X, params)
end

function matrix1d_ideal_HQ_density!(;dmn_eps=1e-6)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_ideal_HQ_density!(A_i, Source, ϕ, tau, X, params)
end

# ── Initial conditions ───────────────────────────────────────────────────────────────────────────
include("profiles_cent1_cent2.jl")

# The dict initializer also accepts the earlier keys: `:nu` for the field `nur`, and the rest-frame
# keys of the earlier second-moment layout, (PiQr | Mr, PiQperp | Mperp, PiQ | Meta):
# piphiphi_Q = PiQperp, pietaeta_Q = -(PiQr + PiQperp), piB_Q = PiQ (`:Meta` meant the trace there).
# A key in the new names wins.
function initialize_fields(init_fun_dict::AbstractDict, field_initializer, grid_params)
    d = Dict{Symbol,Any}(init_fun_dict)
    haskey(d, :nur) || !haskey(d, :nu) || (d[:nur] = d[:nu])
    get2(a, b) = haskey(d, a) ? d[a] : get(d, b, nothing)
    fp = get2(:PiQperp, :Mperp); fr = get2(:PiQr, :Mr); fB = get2(:PiQ, :Meta)
    haskey(d, :piphiphi_Q) || fp === nothing || (d[:piphiphi_Q] = fp)
    haskey(d, :pietaeta_Q) || fr === nothing || fp === nothing || (d[:pietaeta_Q] = x -> -(fr(x) + fp(x)))
    haskey(d, :piB_Q) || fB === nothing || (d[:piB_Q] = fB)
    return invoke(initialize_fields, Tuple{Any,Any,Any}, d, field_initializer, grid_params)
end

# ── Detectors ────────────────────────────────────────────────────────────────────────────────────
include("detector_struct.jl")

export detector_dict, detector_collection, free_hadron, HQdiffusion, conformal_aSB
export oneshoot_debug, spectra_analytic, multiplicity_analytic
