# Charm first moment on a prescribed 1+1D cylindrical background.
#
# The background (T, u) is prescribed, not evolved. `background_fields(tau, r)` returns
#     T, ur, dtT, drT, drur, dtur
# and may depend on both tau and r; "const_BG" in the file name is historical. At c_M = 0 the
# charm does not back-react on the medium, so this is equivalent to co-evolving the two.
#
# Fields: phi = (alpha, nu^r).
# Equations: Nabla.(n u + nu) = 0;  tau_n Delta D nu + nu + kappa Delta grad(alpha) = 0
# (the grad-alpha-only closure). The entries agree with a symbolic covariant derivation.
# The thermodynamically consistent first moment is in HQ_const_BG_consistent.jl.




function matrix1d_visc_HQ_BG!(;dmn_eps=1e-6, background_fields, α_max=200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields, α_max=α_max)
end


# Stable on the Gubser-flow bump test.
# dmn_eps regulates dn/dalpha where n -> 0; both methods default to 1e-6.
function matrix1d_visc_HQ_BG!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-6,background_fields= nothing, α_max=200.0)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    
    thermo = thermodynamic(T,α_safe,params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu+= dmn_eps

    κ = diffusion_hadron(T,α_safe,params.eos,params.diffusion) #diffusion coefficient for hadrons
    # τ_n with the diagnostic knobs (HQ_TAUN_SCALE, FLUIDUM_HQ_MOMENTUM_DIMS) applied, as in every
    # charm matrix. At their defaults this is the bare τ_n.
    tauDiffusion=hq_taun(T,α_safe,params.eos,params.diffusion)

    (At,Ax, source)=one_d_viscous_matrix_fugacity_BG_only(ϕ,tau,X[1],ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,κ)

        
    Ainv= inv(At)


    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

    end 



function one_d_viscous_matrix_fugacity_BG_only(X,tau,r,ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,κ)
    # the charm first moment on a prescribed medium, every row multiplied by u^tau^3 r tau
    At=SMatrix{2,2}(((1  .+ ur .^2) .^2 .*r .*tau .*dn_dmu,

κ .*ur .*(1  .+ ur .^2) .^2 .*r .*tau,

ur .*(1  .+ ur .^2) .*r .*tau,

(1  .+ ur .^2) .^2 .*tauDiffusion .*r .*tau))
        #########################################################################################################################################################################
        
    Ax=SMatrix{2,2}(
        (ur .*(1  .+ ur .^2) .^1.5 .*r .*tau .*dn_dmu,

  κ .*(1  .+ ur .^2) .^2.5  .*r .*tau,

(1  .+ ur .^2) .^1.5 .*r .*tau,

ur .*(1  .+ ur .^2) .^1.5 .*tauDiffusion .*r .*tau)
    )
        #########################################################################################################################################################################
    
    source=SVector{2}(
    (X[2] .*(sqrt.(1  .+ ur .^2) .*tau  .+ ur .^2 .*sqrt.(1  .+ ur .^2) .*tau  .+ r .*(ur  .+ ur .^3  .+ dtur .*tau))  .+ (1  .+ ur .^2) .*n .*(ur .*sqrt.(1  .+ ur .^2) .*tau  .+ r .*(1  .+ ur .^2  .+ dtur .*ur .*tau  .+ drur .*sqrt.(1  .+ ur .^2) .*tau))  .+ (1  .+ ur .^2) .*(dtT  .+ dtT .*ur .^2  .+ drT .*ur .*sqrt.(1  .+ ur .^2)) .*r .*tau .*dn_dT,

(sqrt.(1  .+ ur .^2)  .- ur .*(dtur .*tauDiffusion  .+ dtur .*ur .^2 .*tauDiffusion  .+ ur .*sqrt.(1  .+ ur .^2) .*( .-1  .+ drur .*tauDiffusion))) .*r .*X[2] .*tau)
    )
        
    return (At,Ax, source)
end


