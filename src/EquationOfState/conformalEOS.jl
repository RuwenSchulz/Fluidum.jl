"""
    ConformalEOS(; g_eff = 40.0)

Conformal (Stefan-Boltzmann) medium: `P = a_SB T^4` with `a_SB = g_eff π²/90`, in Fluidum's
units — T in GeV, P in GeV/fm³ (hence the `fmGeV^3`), s = dP/dT in fm⁻³ — the units of
`FluiduMEoS`. Note FluiduMEoS's high-T limit is the 2+1-flavour value g_eff = 47.5
(`95π²/180 = 47.5 π²/90` in `num`); the default 40 here is kept because existing callers rely on it.

Methods (checked 2026-10-13, P3b of Projects/FluidumValidation/MATRICES_PROGRAMME.md, gates in
Projects/FluidumValidation/gate_two_way.jl): `pressure`, `thermodynamic(T, eos)` (P, P', P''),
`thermodynamic_perturbation(T, eos)` (… and P''', for `pressure_derivative(T, Val(3), eos)`) and
`thermodynamic(T, α, eos)` — the medium carries no charge, so this is the α-INDEPENDENT
two-variable form with every α-derivative exactly zero (until 2026-10-13 it had no method, so a
two-variable call MethodError'd). Everything the viscosity models need is P' (entropy) and P''.
"""
Base.@kwdef struct ConformalEOS{T} <: EquationOfState
    g_eff::T = 40.0
end

const π2 = π^2
@inline a_SB(eos::ConformalEOS) = (π2/90) * eos.g_eff

"""
    conformal_aSB(eos::ConformalEOS)

The Stefan-Boltzmann coefficient in Fluidum's units, `P = conformal_aSB(eos) T^4` with P in
GeV/fm³ and T in GeV — i.e. `a_SB(eos) * fmGeV^3`. This is the `aSB` argument of the CONFORMAL
generated two-way functions (`*_conformal_generated` in Matrix/two_way_generated.jl).
"""
@inline conformal_aSB(eos::ConformalEOS) = a_SB(eos) * fmGeV^3

@inline @fastmath pressure(T, eos::ConformalEOS) = a_SB(eos) * T^4 * fmGeV^3

thermodynamic(T::N, eos::ConformalEOS{S}) where {N,S} = begin
    # float(): ConformalEOS(40) with an integer T would otherwise try TT = Int (InexactError)
    TT = float(promote_type(N,S))

    a = TT(a_SB(eos)) * TT(fmGeV^3)
    T2 = T*T
    T3 = T2*T
    T4 = T2*T2

    p    = a * T4
    dp   = 4a * T3
    d2p  = 12a * T2

    Thermodynamic{TT,1,1}(p, (dp,), (d2p,))
end

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
