# ─────────────────────────────────────────────────────────────────────────────────────────────
# ⛔ RETRACTED IN PLACE 2026-10-13 (Projects/FluidumValidation/MATRICES_PROGRAMME.md, P2).
# DEAD: the wrapper builds a NON-FINITE matrix at every state (FLUIDUM_FULL_AUDIT.md §1.2: the
# charm EoS returns an incomplete Hessian, ddP/dT² = ddP/dα² = 0, so A_t is singular), and the
# `invert_n_for_alpha` it calls is off by the fm⁻³ ↔ GeV³ factor 130 (§1.3, gate X4).
# The kernel (two-way: charm pressure in T^{mn}) is UNVERIFIED -- no gate has ever compared it with
# a derivation (its source notebook is Notebooks/ideal/ideal_r_tau_br_p_n.nb). Its mode (:ideal_HQ_br_p_density) is registered in
# MainFluidum/FluidumEvolutionModes.jl but nothing runs it. Kept, not repaired, because this repo
# keeps its wrong turns; the two-way charm system is being derived from scratch (P3a,
# Julia/tools/derive_two_way*.wls). Do NOT "fix" the wrapper by binding the missing names -- see
# FLUIDUM_FULL_AUDIT.md §1.1-1.3 for why the obvious patch is wrong.
# ─────────────────────────────────────────────────────────────────────────────────────────────
#
# The matrices are based on the notebook: ideal_r_tau_br_p_n.nb
# The fields are stored in the order : T, ur, n 
# HQ back reaction and pressure feedback included
#

function matrix1d_ideal_HQ_br_p_density!(A_i, Source, ϕ, t, X, params;free=true)

    #@show t, ϕ
    dP_dT = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dP_dTdT = pressure_derivative(ϕ[1], Val(2), params.eos)



        α = invert_n_for_alpha(ϕ[1], ϕ[3]; m=params.diffusion.mass)

        therm = thermodynamic(ϕ[1], α, params.eos.hadron_list)

     
        dPhq_dT, dPhq_dalpha = therm.pressure_derivative
        ddPhq_dTdT, ddPhq_dTdalpha, ddPhq_dalphadalpha = therm.pressure_hessian



    (At, Ax, source) = one_d_ideal_matrix_ruwen6(ϕ, t, X[1],dP_dT, dP_dTdT,dPhq_dT,dPhq_dalpha, ddPhq_dTdT, ddPhq_dTdalpha,ddPhq_dalphadalpha)

    Ainv = inv(At)

    A_mul_B!(A_i[1], Ainv, Ax)

    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen6(X, tau, r, dP_dT, dP_dTdT,dPhq_dT,dPhq_dalpha, ddPhq_dTdT, ddPhq_dTdalpha,ddPhq_dalphadalpha)

    At = SMatrix{3,3}((( .-dPhq_dalpha .^2  .+ X[1] .*(2  .+ X[2] .^2) .*dPhq_dalpha .*ddPhq_dTdalpha  .+ X[1] .*(X[2] .^2 .*dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT .*ddPhq_dalphadalpha  .+ X[2] .^2 .*ddPhq_dalphadalpha .*dPhq_dT  .- X[1] .*ddPhq_dTdalpha .^2  .- X[1] .*X[2] .^2 .*ddPhq_dTdalpha .^2  .+ X[1] .*ddPhq_dalphadalpha .*ddPhq_dTdT  .+ X[1] .*X[2] .^2 .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./(X[1] .*ddPhq_dalphadalpha),

(X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*dP_dTdT .*ddPhq_dalphadalpha  .+ ddPhq_dalphadalpha .*dPhq_dT  .+ dPhq_dalpha .*ddPhq_dTdalpha  .- X[1] .*ddPhq_dTdalpha .^2  .+ X[1] .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./ddPhq_dalphadalpha,

0,

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

(X[2] .*X[3]) ./sqrt.(1  .+ X[2] .^2),

(X[1] .^2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdalpha) ./ddPhq_dalphadalpha,

(X[1] .*(dPhq_dalpha  .+ X[1] .*X[2] .^2 .*ddPhq_dTdalpha)) ./ddPhq_dalphadalpha,

sqrt.(1  .+ X[2] .^2)))
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}(((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*dP_dTdT .*ddPhq_dalphadalpha  .+ ddPhq_dalphadalpha .*dPhq_dT  .+ dPhq_dalpha .*ddPhq_dTdalpha  .- X[1] .*ddPhq_dTdalpha .^2  .+ X[1] .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./ddPhq_dalphadalpha,

(dPhq_dalpha .^2  .+ X[1] .*( .-1  .+ X[2] .^2) .*dPhq_dalpha .*ddPhq_dTdalpha  .+ X[1] .*((1  .+ X[2] .^2) .*dP_dT .*ddPhq_dalphadalpha  .+ X[1] .*X[2] .^2 .*dP_dTdT .*ddPhq_dalphadalpha  .+ ddPhq_dalphadalpha .*dPhq_dT  .+ X[2] .^2 .*ddPhq_dalphadalpha .*dPhq_dT  .- X[1] .*X[2] .^2 .*ddPhq_dTdalpha .^2  .+ X[1] .*X[2] .^2 .*ddPhq_dalphadalpha .*ddPhq_dTdT)) ./(X[1] .*ddPhq_dalphadalpha),

0,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*(dP_dT  .+ dPhq_dT)) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*(dP_dT  .+ dPhq_dT),

X[3],

(X[1] .*( .-dPhq_dalpha  .+ X[1] .*(1  .+ X[2] .^2) .*ddPhq_dTdalpha)) ./ddPhq_dalphadalpha,

(X[1] .^2 .*X[2] .*sqrt.(1  .+ X[2] .^2) .*ddPhq_dTdalpha) ./ddPhq_dalphadalpha,

X[2]))
    #########################################################################################################################################################################

    source = SVector{3}(((X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

(X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*(dP_dT  .+ dPhq_dT)) ./(r .*tau),

X[3] .*(X[2] ./r  .+ sqrt.(1  .+ X[2] .^2) ./tau)))
         return (At,Ax, source)
end

