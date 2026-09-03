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
    tauq = 0
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
    norm = 0
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
    norm = 0
    for i in x.hadron_list.particle_list
        m = i.Mass
        q = i.Nc + i.Nac
        n = free_hadron(T,μ,1,q; m = m)[1]
        norm += q^2*n
    end
    return norm;
end


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



