Base.@kwdef struct Heavy_Quark{T,S} <:EquationOfState 
    mass::T = 1.5
    eos::FluiduMEoS{T} = FluiduMEoS()
    hadron_list::HadronResonaceGas_ccbar{S,T}  
end

function Heavy_Quark(particle_list, ccbar)
    list = HadronResonaceGas_ccbar(particle_list, ccbar)
    return Heavy_Quark(hadron_list = list)
end

@inline @fastmath pressure(T, x::Heavy_Quark) = pressure(T, x.eos)
@inline pressure_derivative(T, v::Val{1}, x::Heavy_Quark) = pressure_derivative(T, v, x.eos)
@inline pressure_derivative(T, v::Val{2}, x::Heavy_Quark) = pressure_derivative(T, v, x.eos)

function free_hadron(T,μ,deg,q; m = 1.5)
    b1 = besselk1x(m/T)
    b2 = besselkx(2,m/T) 
    b3 = b1+4/(m/T)*b2

    ex = exp(q*μ - m/T)  
    density = deg*(T /(2 *π^2)*m^2* ex* b2)* fmGeV3; #fm-3
    densityDerT = deg*((ex*m^2*(m*b1 +2*T*b2 + m*b3))/(4*π^2*T))* fmGeV3 ; #(*fm^-3/GeV)
    
    densityDerμ= density/T+0.0001;
    densityDerTDerμ =deg*(ex*m^2*(m*b1 + 2*T*b2 + m*b3)/(4*π^2*T))* fmGeV3 ;
    return (density,densityDerT,densityDerμ,densityDerTDerμ)
end

# struct Const_DsT{T} <:DiffusionModel
#     DsT::T
# end

# struct Linear_DsT{T} <:DiffusionModel
#     slope::T
#     offset::T
# end

# struct HQdiffusion{L,T}<:Diffusion
#     DsT_model::L
#     mass::T
# end

# HQdiffusion(DsT::M, mass::M) where {M} = HQdiffusion{Const_DsT{M},M}(Const_DsT(DsT),mass)
# HQdiffusion(slope::M, offset::M, mass::M) where {M} = HQdiffusion{Linear_DsT{M},M}(Linear_DsT(slope, offset), mass)


# DsT(m::Const_DsT,T) = m.DsT
# function DsT(m::Linear_DsT,T;Tfo = 0.156)
#     if T > Tfo
#         return m.slope*T + m.offset
#     else
#         return m.slope*Tfo + m.offset
#     end
# end 

# DsT(y::HQdiffusion, T) = DsT(y.DsT_model, T)
#DsT(y::HQdiffusion, T) = DsT(y.DsT_model, T)

@inline function DsT(m::ConstDiffusion,T) 
    m.DsT
end

@inline function DsT(m::LinearDiffusion,T;Tfo = 0.156)
    if T > Tfo
        return m.slope*T + m.offset
    else
        return m.slope*Tfo + m.offset
    end
end 


# 2026-09-02 FIX: the body referenced an undefined `x` (`x.DsT(y,T)`), so ANY call with a
# ConstDiffusion / LinearDiffusion threw `UndefVarError: x not defined in Fluidum`. This is the
# error `Projects/NumericBackendAudit/AUDIT_REPORT.md` recorded against `oneshoot_debug`: the
# causality checker calls `diffusion(T,n,params.diffusion)` and `τ_diffusion(T,params.diffusion)`,
# and both were dead. `DsT` is a plain function, never a field of the coefficient struct — cf.
# `diffusion_hadron` just below, the same expression with `normalization(T,μ,x)` in place of `n`.
function diffusion(T,n,y::Diffusion)
    density = n #fm-^3
    κ = DsT(y,T)/T*density/fmGeV #fm^-2
end

function diffusion_hadron(T,μ,x::Heavy_Quark,y::Diffusion)
    κ = DsT(y,T)/T*normalization(T,μ,x)/fmGeV
end


# 2026-09-02 FIX: the same undefined-`x` defect (`x.mass`, `x.DsT(y,T)`). The mass lives on the
# diffusion coefficient itself (ConstDiffusion.mass / LinearDiffusion.mass). With `y` substituted
# this method is ALGEBRAICALLY IDENTICAL to `τ_diffusion_hadron` on a single unit-charge species of
# mass `y.mass`: both reduce to  τ_n = D_sT·m³(2b₁-3b₃+b₅)/(48 T⁴ b₂)/fmGeV, the degeneracy
# cancelling against `normalization`. Gated as D1 in Projects/FluidumValidation/bench_fluidum_attractor.jl.
function τ_diffusion(T,y::Diffusion) 
    m = y.mass
  
    b2 = besselkx(2,m/T)*exp(-m/T)
    b1 = besselkx(1,m/T)*exp(-m/T)
   
    b3 = b1+4/(m/T)*b2  
    b4 = b2 + 6/(m/T)*b3  
    b5 = b3+8/(m/T)*b4

  return ((2*π*DsT(y,T)) *m^3/T^3/(96*π*T)*(2*b1 - 3*b3 +b5)/b2)/fmGeV; #
  
end


function τ_diffusion_hadron(T,μ,x::Heavy_Quark,y::Diffusion)
    tauq = zero(float(promote_type(typeof(T), typeof(μ))))
    for i in x.hadron_list.particle_list
        m = i.Mass
        q = i.Nc + i.Nac
        deg = i.Degeneracy
        b2 = besselkx(2,m/T)*exp(-m/T)
        b1 = besselk1x(m/T)*exp(-m/T)
        b3 = b1+4/(m/T)*b2
        b4 = b2 + 6/(m/T)*b3
        b5 = b3+8/(m/T)*b4
        ex = exp(q*μ)
        # 2026-07-21 FIX: the numerator carries `deg` so it cancels the `deg` in `normalization`
        # (= Σ q²·free_hadron ∝ deg). τ_n is a ratio of moments of the SAME distribution ⇒ the
        # degeneracy MUST cancel — this is Eq. (30) of Capellino et al. 2205.07692, τ_n = D_s·I₃₁/(T·P₀),
        # and I₃₁, P₀ both ∝ g. WITHOUT `deg` here the single-species charm quark (g=6) came out
        # τ_n = τ_n^bare/6 (superluminal above T≈0.48; contradicts the reference, the FP derivation,
        # LangevInMedium.tau_n_main3, and Fluidum's own first-moment τ_diffusion). Harmless for the
        # original pseudoscalar-D-meson list (deg=1); a genuine 6× error only for the charm quark.
        tauq += (deg*(2*π*DsT(y,T))/(192*π^3*T^3)*q^2*m^5*ex*(2*b1 - 3*b3 +b5));
    end
  return tauq/normalization(T,μ,x)*(fmGeV^2); #
end


# ── THE charm relaxation times, with the two diagnostic knobs applied ONCE (2026-10-13) ──────────
# Until this date each matrix file built τ_n / τ_M / η_M / h itself and the knobs disagreed:
# HQ_TAUN_SCALE was ignored by the 1-D first moment; FLUIDUM_HQ_MOMENTUM_DIMS=2 rescaled τ_M in the
# (now retracted) shipped second moment but never in the consistent one; the 2-D m2 driver paired an
# UNrescaled τ_n with an h built from the rescaled one (ratio 0.9156 at DIMS=2, breaking τ_n = D_s h/T).
# Every charm matrix now takes them from here. At the production defaults (scale 1.0, DIMS 3) the
# values are bit-identical to what each file computed before.
#
# FLUIDUM_HQ_MOMENTUM_DIMS=2 (2026-08-21): match the charm sector in the TWO-dimensional momentum
# measure. The Langevin ensembles of LangevinPaper1 evolve transverse momenta only, so their
# equilibrium is the 2-D Jüttner p dp e^{-E/T} and their ℓ-th moment relaxes at λ_ℓ(2D)·η_D, not at
# the 3-D λ_ℓ = K_{ℓ+1}/K_{ℓ+2} that τ_n = D_s z K₃/K₂ and τ_M = D_s z K₄/(2K₃) encode. In that
# measure every moment ∫_M^∞ E^n e^{-E/T} dE is an incomplete Gamma of integer order, so the ratios
# are closed forms (checked against quadrature to 1e-9 at z = 2…20):
#     λ₁(2D) = z(z+1)/(z²+3z+3),      λ₂(2D) = z(z²+3z+3)/(z³+6z²+15z+15),
# and τ_n^{2D} = τ_n·(K₂/K₃)/λ₁(2D)  (5–12 % shorter over z = 3.5–10),  τ_M^{2D} = τ_M·(K₃/K₄)/λ₂(2D).
# η_M stays tied to τ_n (η_M = T τ_n/2); κ = D_s n is dimension-independent. The 5/3 shear/bulk
# tie and the tensor structure of the matrices are NOT switched (they are the 3-spatial-dimension
# projection the solvers are written in). DIAGNOSTIC knob.

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

The charm transport coefficients of the first and second moment, tied to ONE τ_n:

    τ_M = τ_n K₄K₂/(2K₃²)   (use_NR_tauM: τ_n/2),    η_M = T τ_n/2,
    D_s = DsT/T (fm),        h = τ_n T/D_s           (so τ_n = D_s h/T exactly),

with `HQ_TAUN_SCALE` and `FLUIDUM_HQ_MOMENTUM_DIMS` applied (z = m/T, m = `diffusion.mass`). The
FP→hydro matching (Tex/HydroFPderivation/FP_Hydro_matching.tex) defines τ_n and τ_M in the same
convention, so the RATIO τ_M/τ_n = ½K₄K₂/K₃² (≈ 0.550 at z = 9.615) is what is physical; tying τ_M
to τ_n keeps that ratio whatever τ_n's normalisation.
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


# ── THE CONSISTENT SECOND MOMENT'S EOS COEFFICIENTS FOR A CHARM LIST (2026-10-13) ─────────────────
# The consistent second moment (Tex/LangevinPaper1/M2_CONSISTENT_DERIVATION.md) needs, besides h, four
# EOS objects of the charm gas: I31 (its T-log-derivative A/B multiplies D ln T in the trace drive),
# A = I41/P0 (the a·ν completion terms), C = I52/I42 (τ_M = D_s C/(2T) and D ln C) — derived for ONE
# Boltzmann species. For a LIST of species i (mass m_i, degeneracy g_i, charm charge q_i) at a common
# (T, α) with ONE ν and ONE δM field the moments are species sums, with the weights the existing
# aggregates already carry (n = Σ q n_i; ν channel δf_i ∝ q_i f0_i ⇒ h = Σ q² n_i h_i / Σ q² n_i, the
# weighting for which D_s h/T IS τ_diffusion_hadron; π channel δf_i ∝ q_i f0_i likewise):
#     I31 = n T h          ⇒ A/B = ∂ln I31/∂ln T|_α = T ∂_T n/n + 1 + T ∂_T h/h,   ∂_α ln I31 = (Σ q² n_i)/n
#     A   = Σ q² n_i (m_i² + 5T h_i) / Σ q² n_i  =  ⟨m²⟩ + 5T h,         ⟨m²⟩ = Σ q² n_i m_i² / Σ q² n_i
#     C   = Σ q² I42_i C_i / Σ q² I42_i,   I42_i = T² n_i h_i,   C_i = m_i K4(z_i)/K3(z_i)
# with every T-derivative analytic at fixed α. For one species of mass m all of these are the
# single-quark closed forms (A/B = 5 + zK2/K3, A = m² + 5Th, C = m K4/K3).
#
# 🔴 WHY THIS EXISTS. Until 2026-10-13 the m2 DRIVERS (HQ_const_BG_consistent_m2.jl, HQ_2p1d_BG_m2.jl)
# took A/B, d ln C/d ln T, A and C from the SINGLE-QUARK forms at m = diffusion.mass, and the two-way
# WRAPPERS (two_way.jl) took A/B from the list but C = h K4K2/K3² and A = m² + 5Th at diffusion.mass.
# Against a direct momentum quadrature over the 77-state production list (T = 0.156–0.45 GeV):
# drivers A/B −22…−26 %, d ln C/d ln T −24…+26 %; both C (hence τ_M) +1.6…+6.6 %, A −27…−37 %;
# wrappers' d ln C/d ln T +5…+12 % (MATRICES_PROGRAMME.md Log 2026-10-13). Identical for one species.
#
# ⚠ UNIT-CHARGE ASSUMPTION (the derivation's). The one q = 2 state of the production list (J/ψ) enters
# with the weights above; a charge-correct second moment would also carry ∂_α C, ∂_α h, and Σ q I31_i
# instead of nTh. Size at production list loadings: MATRICES_PROGRAMME.md Log 2026-10-13.

"""
    with_m2_mixture(f, on::Bool)

Run `f()` with `HQ_M2_MIXTURE[] = on` and restore it afterwards (closures read the Ref when they are
BUILT, so build them inside `f`). `on = false` reproduces the pre-2026-10-13 list coefficients.
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

True when the charm `list` is ONE unit-charge species of mass `mq` (LP1's and O+O's `charm.data`):
the second-moment coefficients are then the single-quark closed forms and every m2 matrix keeps
its pre-2026-10-13 code path bit for bit.
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
(keeping the diagnostic knobs' factors). The one implementation behind the m2 drivers AND the two-way
wrappers. Zero allocations.
"""
function hq_m2_mixture(T, α, list::HadronResonaceGas_ccbar, mq)
    R = float(promote_type(typeof(T), typeof(α)))
    # common scale e^{-sh}, sh = max_i (q_i α − z_i): every weight is a RATIO, so nothing underflows
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


struct QGPViscosity{T}<:ShearViscosity
    ηs::T
    Cs::T
end


@inline function bulk_viscosity(T,entropy,y::SimpleBulkViscosity{N}) where {N}
    y.ζs/(1+((T-0.175)/0.024)^2)*invfmGeV*entropy
 end

@inline function τ_bulk(T,entropy,dtdtp,y::SimpleBulkViscosity{N}) where {N}
     cs2= entropy/(T*dtdtp)
     bulk_viscosity(T,entropy,y)/(T*entropy*y.Cζ)*1/(1/3-cs2)^2+0.1
 end


 @inline function bulk_viscosity(T,y::ZeroBulkViscosity)  
    zero(promote_type(typeof(T)))
end

@inline function bulk_viscosity(T,entropy,y::ZeroBulkViscosity)  
    zero(promote_type(typeof(T)))
end

@inline function τ_bulk(T,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end

@inline function τ_bulk(T,entropy,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end

@inline function τ_bulk(T,entropy,dtdtp,y::ZeroBulkViscosity)  
    one(promote_type(typeof(T)))
end 


@inline function viscosity(T,entropy,y::QGPViscosity{N}) where {N}
    y.ηs*entropy*invfmGeV
 end

# 2026-09-02: the two shear models carried DISJOINT method sets — SimpleShearViscosity only took a
# `Thermodynamic`, QGPViscosity only a scalar entropy — so `matrix2d_visc!` (which passes the whole
# Thermodynamic) MethodError'd with QGPViscosity, the model every 1-D production run uses, while
# `matrix1d_visc_HQ!` and the new 2+1D solvers (which pass the scalar) MethodError'd with
# SimpleShearViscosity, the one Fluidum's own `2d viscous` testset uses.  Each entry point worked
# with exactly the model the other could not take.  The formulas are IDENTICAL — η = ηs·s/ħc and
# τ_s = η/(T s C_s) = ηs/(ħc T C_s) — so the fix is the four missing dispatches, not new physics.
@inline function viscosity(T,entropy::Number,y::SimpleShearViscosity{N}) where {N}
    y.ηs*entropy*invfmGeV
end
@inline function τ_shear(T,entropy::Number,y::SimpleShearViscosity{N}) where {N}
    y.ηs*invfmGeV/(T*y.Cs)
end
@inline function viscosity(T,x::Thermodynamic{N,1,1},y::QGPViscosity{N}) where {N}
    y.ηs*(@inbounds x.pressure_derivative[1])*invfmGeV
end
@inline function τ_shear(T,x::Thermodynamic{N,1,1},y::QGPViscosity{N}) where {N}
    y.ηs*invfmGeV/(T*y.Cs)
end


 
@inline function τ_shear(T,entropy,y::QGPViscosity{N}) where {N}
    viscosity(T,entropy,y)/((T*entropy)*y.Cs)
end


@inline function viscosity(T,entropy,y::ZeroViscosity) 
    zero(promote_type(typeof(T),typeof(entropy)))
 end


 
@inline function τ_shear(T,entropy,y::ZeroViscosity) 
   one(promote_type(typeof(T),typeof(entropy)))
end



function normalization(T,μ,x::Heavy_Quark)
    norm = zero(float(promote_type(typeof(T), typeof(μ))))
    for i in x.hadron_list.particle_list
        m = i.Mass
        deg = i.Degeneracy
        q = i.Nc + i.Nac
        n = free_hadron(T,μ,deg,q; m = m)[1]
        norm += q^2*n
    end 
    return norm;    
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


# The ÷g_hq that `HQ_const_BG_2nd_moment.jl` once divided τ_M/η_M by (removed 2026-07-21). No
# matrix calls it any more; kept because Julia/tools/rhs_compare.jl still reads it to reproduce the
# old convention.
second_moment_transport_normalization(T,μ,x) = one(promote_type(typeof(T), typeof(μ)))

function second_moment_transport_normalization(T,μ,x::Heavy_Quark)
    bare_norm = normalization_without_degeneracy(T,μ,x)
    bare_norm == 0 && return one(promote_type(typeof(T), typeof(μ)))
    return normalization(T,μ,x) / bare_norm
end




function τ_diffusion(T,x::ZeroDiffusion) 
    return 1.0
end

function τ_diffusion(T,n, x::ZeroDiffusion) 
    return 1.0
end

function τ_diffusion_hadron(T,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 1.0
end

function τ_diffusion_hadron(T,n,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 1.0
end


function diffusion(T,n,x::ZeroDiffusion) 
    return 0.0
end

function diffusion(T,x::ZeroDiffusion) 
    return 0.0
end


function diffusion_hadron(T,n,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 0.0
end

function diffusion_hadron(T,μ,x::Heavy_Quark,y::ZeroDiffusion) 
    return 0.0
end



