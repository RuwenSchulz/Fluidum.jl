# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, mu, nur

function matrix1d_visc_HQ_BG_second_moment!(;dmn_eps=1e-6, background_fields)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_second_moment!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
function matrix1d_visc_HQ_BG_second_moment!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-10,background_fields= nothing)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])
    
    thermo = thermodynamic(T,ϕ[1],params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu+= dmn_eps

    κ = diffusion_hadron(T,ϕ[1],params.eos,params.diffusion) #diffusion coefficient for hadrons
    tauDiffusion=τ_diffusion_hadron(T,ϕ[1],params.eos,params.diffusion) #tau diffusion for hadrons
    tau_phi = tauDiffusion/2
    etaPhi = 1.5 * T * tau_phi
    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )
    #(At,Ax, source)=one_d_viscous_HQ_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix_fugacity_BG_only_second_moment(ϕ,tau,X[1],ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,κ,tau_phi,etaPhi)

    Ainv= inv(At)

    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

    end 



function one_d_viscous_matrix_fugacity_BG_only_second_moment(X,tau,r,ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,κ,tauPhi,etaPhi)
    #these are the matrices from 1d viscous hydro
    At=SMatrix{3,3}(
((1  .+ ur .^2) .^2 .*r .*tau .*dn_dmu,

ur .*(1  .+ ur .^2) .^2 .*κ.*r .*tau,

0,

ur .*(1  .+ ur .^2) .*r .*tau,

(1  .+ ur .^2) .^2 .*tauDiffusion .*r .*tau,

(4 .*etaPhi .*ur .*(1  .+ ur .^2) .^2 .*r .*tau) ./3.,

0,

0,

(tauPhi .*(1  .+ ur .^2) .^3 .*(3  .+ 2 .*ur .^2) .*r .*tau) ./3.)
    )
        #########################################################################################################################################################################
        
    Ax=SMatrix{3,3}(
(ur .*(1  .+ ur .^2) .^1.5 .*r .*tau .*dn_dmu,

(1  .+ ur .^2) .^2.5 .*κ.*r .*tau,

0,

(1  .+ ur .^2) .^1.5 .*r .*tau,

ur .*(1  .+ ur .^2) .^1.5 .*tauDiffusion .*r .*tau,

(4 .*etaPhi .*(1  .+ ur .^2) .^2.5 .*r .*tau) ./3.,

0,

(1  .+ ur .^2) .^2.5 .*r .*tau,

(tauPhi .*ur .*(1  .+ ur .^2) .^2.5 .*(3  .+ 2 .*ur .^2) .*r .*tau) ./3.)
    )
        #########################################################################################################################################################################
    
    source=SVector{3}(
(X[2] .*(sqrt.(1  .+ ur .^2) .*tau  .+ ur .^2 .*sqrt.(1  .+ ur .^2) .*tau  .+ r .*(ur  .+ ur .^3  .+ dtur .*tau))  .+ (1  .+ ur .^2) .*n .*(ur .*sqrt.(1  .+ ur .^2) .*tau  .+ r .*(1  .+ ur .^2  .+ dtur .*ur .*tau  .+ drur .*sqrt.(1  .+ ur .^2) .*tau))  .+ (1  .+ ur .^2) .*(dtT  .+ dtT .*ur .^2  .+ drT .*ur .*sqrt.(1  .+ ur .^2)) .*r .*tau .*dn_dT,

r .*(2 .*ur .*(1  .+ ur .^2) .^2 .*X[3]  .+ (sqrt.(1  .+ ur .^2)  .- ur .*(dtur .*tauDiffusion  .+ dtur .*ur .^2 .*tauDiffusion  .+ ur .*sqrt.(1  .+ ur .^2) .*( .-1  .+ drur .*tauDiffusion))) .*X[2] .*tau),

 .-0.3333333333333333 .*((1  .+ ur .^2) .*(4 .*etaPhi .*ur .*(dtur .*ur  .+ drur .*sqrt.(1  .+ ur .^2)) .*r .*X[2] .*tau  .- 3 .*sqrt.(1  .+ ur .^2) .*r .*X[3] .*tau  .+ 2 .*tauPhi .*(1  .+ ur .^2) .*X[3] .*(2 .*(1  .+ ur .^2) .*r  .- ur .*sqrt.(1  .+ ur .^2) .*tau))))
 )
        
    return (At,Ax, source)
end
