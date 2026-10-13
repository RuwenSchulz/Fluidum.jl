#
# The matrices are based on the notebook: ideal_r_tau_n.nb
# The current is on the EOS of the fluid, i.e. no back reaction of the current, no back feeding of the HQ pressure
# The fields are stored in the order : T, ur, n 
#
# Equations: ideal Nabla_m T^{mn} = 0 and Nabla.(n u) = 0 (no diffusion current). Row 3 has d_tau T
# and d_tau u^r eliminated with rows 1-2. GATED 2026-10-13: M-I3n of Julia/tools/derive_1d_matrices.wls
# (solved form, exact rational arithmetic) and J-In of Projects/FluidumValidation/gate_matrices_1d.jl.
# Mode :ideal_HQ_density; nothing in production runs it.
#
# `dmn_eps` is ACCEPTED AND IGNORED: the mode registry calls every closure as `(; dmn_eps = ...)`,
# and this matrix has no dn/dalpha to regularise (n is the field). The unused `free` keyword the
# direct method carried was removed 2026-10-13 (no caller passed it).

function matrix1d_ideal_HQ_density!(;dmn_eps=1e-6)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_ideal_HQ_density!(A_i, Source, ϕ, tau, X, params)
end


function matrix1d_ideal_HQ_density!(A_i, Source, ϕ, t, X, params)

    #@show t, ϕ
    dP_dT = pressure_derivative(ϕ[1], Val(1), params.eos) #entropy
    dP_dTdT = pressure_derivative(ϕ[1], Val(2), params.eos)


    (At, Ax, source) = one_d_ideal_matrix_ruwen4(ϕ, t, X[1],dP_dT, dP_dTdT)

    Ainv = inv(At)

    A_mul_B!(A_i[1], Ainv, Ax)

    jgemvavx!(Source, Ainv, source)



end


function one_d_ideal_matrix_ruwen4(X, tau, r, dP_dT, dP_dTdT)

    At = SMatrix{3,3}((X[2] .^2 .*dP_dT  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT,

X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

0,

2 .*X[1] .*X[2] .*dP_dT,

(X[1] .*(1  .+ 2 .*X[2] .^2) .*dP_dT) ./sqrt.(1  .+ X[2] .^2),

0,

0,

0,

sqrt.(1  .+ X[2] .^2)))
    #########################################################################################################################################################################

    Ax = SMatrix{3,3}((X[2] .*sqrt.(1  .+ X[2] .^2) .*(dP_dT  .+ X[1] .*dP_dTdT),

(1  .+ X[2] .^2) .*dP_dT  .+ X[1] .*X[2] .^2 .*dP_dTdT,

 .-((X[2] .*X[3] .*dP_dTdT) ./( .-(X[2] .^2 .*dP_dT)  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT)),

(X[1] .*(1  .+ 2 .*X[2] .^2) .*dP_dT) ./sqrt.(1  .+ X[2] .^2),

2 .*X[1] .*X[2] .*dP_dT,

(X[1] .*X[3] .*dP_dTdT) ./( .-(X[2] .^2 .*dP_dT)  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT),

0,

0,

X[2]))
    #########################################################################################################################################################################

    source = SVector{3}(((X[1] .*(r  .+ r .*X[2] .^2  .+ X[2] .*sqrt.(1  .+ X[2] .^2) .*tau) .*dP_dT) ./(r .*tau),

(X[1] .*X[2] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*dP_dT) ./(r .*tau),

(X[1] .*(1  .+ X[2] .^2) .*X[3] .*(r .*sqrt.(1  .+ X[2] .^2)  .+ X[2] .*tau) .*dP_dTdT) ./(r .*tau .*( .-(X[2] .^2 .*dP_dT)  .+ X[1] .*(1  .+ X[2] .^2) .*dP_dTdT))))
         return (At,Ax, source)
end

