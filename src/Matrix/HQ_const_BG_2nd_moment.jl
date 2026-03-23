function matrix1d_visc_HQ_BG_second_moment!(;dmn_eps=1e-6, background_fields)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_second_moment!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
function matrix1d_visc_HQ_BG_second_moment!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-6,background_fields= nothing)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])
    
    thermo = thermodynamic(T,ϕ[1],params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dalpha = thermo.pressure_derivative
    dn_dalpha+= dmn_eps

    kappa = diffusion_hadron(T,ϕ[1],params.eos,params.diffusion) #diffusion coefficient for hadrons
    tauDiffusion=τ_diffusion_hadron(T,ϕ[1],params.eos,params.diffusion) #tau diffusion for hadrons
    mq = 1.5
    #tau_phi = tauDiffusion/2

    tauM = tauDiffusion/2
    etaM = T * tauM

    #etaPhi  = T * tau_phi
    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )
    #(At,Ax, source)=one_d_viscous_HQ_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix_fugacity_BG_only_second_moment(ϕ,tau,X[1],ur,T,dtT,drT,drur,dtur,n,dn_dalpha,dn_dT,tauDiffusion,kappa,tauM,etaM,mq)
    #println("At = ", At)
    #println("Ax = ", Ax)
    #println("source = ", source)

    Ainv= inv(At)

    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

    end 



    function one_d_viscous_matrix_fugacity_BG_only_second_moment(X,tau,r,ur,T,dtT,drT,drur,dtur,n,dn_dalpha,dn_dT,tauDiffusion,kappa,tauM,etaM,mq)
    #these are the matrices from 1d viscous hydro
    At=SMatrix{5,5}(
(sqrt.(1  .+ ur .^2) .*dn_dalpha,

kappa .*ur .*sqrt.(1  .+ ur .^2),

0,

0,

0,

ur ./sqrt.(1  .+ ur .^2),

sqrt.(1  .+ ur .^2) .*tauDiffusion,

(4 .*etaM .*ur .*sqrt.(1  .+ ur .^2)) ./3.,

( .-2 .*etaM .*ur .*r .^2) ./(3. .*sqrt.(1  .+ ur .^2)),

( .-2 .*etaM .*ur .*tau .^2) ./(3. .*sqrt.(1  .+ ur .^2)),

0,

(ur .*sqrt.(1  .+ ur .^2) .*tauDiffusion) ./mq .^2,

(1  .+ ur .^2) .^1.5 .*tauM,

0,

0,

0,

0,

0,

sqrt.(1  .+ ur .^2) .*tauM .*r .^2,

0,

0,

0,

0,

0,

sqrt.(1  .+ ur .^2) .*tauM .*tau .^2)
    )
        #########################################################################################################################################################################
        
    Ax=SMatrix{5,5}(
(ur .*dn_dalpha,

kappa .*(1  .+ ur .^2),

0,

0,

0,

1,

ur .*tauDiffusion,

(4 .*etaM .*(1  .+ ur .^2)) ./3.,

( .-2 .*etaM .*r .^2) ./3.,

( .-2 .*etaM .*tau .^2) ./3.,

0,

((1  .+ ur .^2) .*tauDiffusion) ./mq .^2,

ur .*(1  .+ ur .^2) .*tauM,

0,

0,

0,

0,

0,

ur .*tauM .*r .^2,

0,

0,

0,

0,

0,

ur .*tauM .*tau .^2)
    )
        #########################################################################################################################################################################
    
    source=SVector{5}(
(n .*(drur  .+ (dtur .*ur) ./sqrt.(1  .+ ur .^2)  .+ ur ./r  .+ sqrt.(1  .+ ur .^2) ./tau)  .+ X[2] .*(1 ./r  .+ (ur  .+ ur .^3  .+ dtur .*tau) ./((1  .+ ur .^2) .^1.5 .*tau))  .+ (drT .*ur  .+ dtT .*sqrt.(1  .+ ur .^2)) .*dn_dT,

(1  .- (dtur .*ur .*tauDiffusion) ./sqrt.(1  .+ ur .^2)  .+ drur .*( .-1  .+ 1 ./(1  .+ ur .^2)) .*tauDiffusion) .*X[2]  .+ ((drur .*ur  .+ dtur .*sqrt.(1  .+ ur .^2)) .*tauDiffusion .*X[3]) ./mq .^2,

(4 .*etaM .*ur .*( .-drur  .- (dtur .*ur) ./sqrt.(1  .+ ur .^2)) .*X[2]) ./3.  .+ (1  .+ ur .^2) .*( .-(n .*T)  .+ X[3]),

(r .*( .-3 .*n .*T .*r  .+ (2 .*etaM .*ur .*(dtur .*ur  .+ drur .*sqrt.(1  .+ ur .^2)) .*r .*X[2]) ./(1  .+ ur .^2) .^1.5  .+ 3 .*( .-2 .*ur .*tauM  .+ r) .*X[4])) ./3.,

(tau .*(2 .*etaM .*ur .*(dtur .*ur  .+ drur .*sqrt.(1  .+ ur .^2)) .*X[2] .*tau  .- 3 .*(1  .+ ur .^2) .*(2 .*(1  .+ ur .^2) .*tauM .*X[5]  .+ n .*T .*sqrt.(1  .+ ur .^2) .*tau  .- sqrt.(1  .+ ur .^2) .*X[5] .*tau))) ./(3. .*(1  .+ ur .^2) .^1.5))
    )
    return (At,Ax, source)
end

