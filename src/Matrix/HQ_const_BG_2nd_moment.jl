function matrix1d_visc_HQ_BG_second_moment!(;dmn_eps=1e-6, background_fields, α_max=200.0, use_NR_tauM::Bool=false)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_second_moment!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields, α_max=α_max, use_NR_tauM=use_NR_tauM)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
function matrix1d_visc_HQ_BG_second_moment!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-6,background_fields= nothing, α_max=200.0, use_NR_tauM::Bool=false)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    
    thermo = thermodynamic(T,α_safe,params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dalpha = thermo.pressure_derivative
    dn_dalpha+= dmn_eps
    #dn_dT += dmn_eps
    kappa = diffusion_hadron(T,α_safe,params.eos,params.diffusion) #diffusion coefficient for hadrons
    Ds = DsT(params.diffusion, T) / T / fmGeV
    mq = params.diffusion.mass

    taun=τ_diffusion_hadron(T,α_safe,params.eos,params.diffusion) #tau diffusion for hadrons
    #@show tau
    z = mq / T
    transport_norm = second_moment_transport_normalization(T, α_safe, params.eos)

    tauM_rel = Ds / (2) * (6 *z *besselk(1, z) + (z^2 + 24) * besselk(2, z))/(z*besselk(1,z) + 4 *besselk(2, z))
    etaM_rel = (Ds*T/2) * (4 + z *besselk(1, z)/besselk(2, z))

    tauM = use_NR_tauM ? taun / 2 : tauM_rel / transport_norm
    etaM = use_NR_tauM ? T * taun / 2 : etaM_rel / transport_norm

    cM   = 0 #Ds / T                # c_M  = D_s/T    (exact, no NR change)

    #(At,Ax, source)=one_d_viscous_HQ_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix_fugacity_BG_only_second_moment(ϕ,tau,X[1],ur,T,dtT,drT,drur,dtur,n,dn_dalpha,dn_dT,taun,kappa,tauM,etaM,mq,cM)

    Ainv= inv(At)

    if !all(isfinite, Ainv)
        @warn "Non-finite inv(At) detected" tau r=X[1] T ur n dn_dalpha taun kappa tauM etaM cM det_At=det(At) maxlog=5
    end

    A_mul_B!(A_i[1], Ainv,Ax)

    jgemvavx!(Source, Ainv,source)

    end 

    function one_d_viscous_matrix_fugacity_BG_only_second_moment(X,tau,r,ur,T,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,kappa,tauM,etaM,mq,cM)
    #these are the matrices from 1d viscous hydro
    At=SMatrix{5,5}(
(sqrt.(1  .+ ur .^2) .*dn_dmu,

kappa .*ur .*sqrt.(1  .+ ur .^2),

0,

0,

0,

ur ./sqrt.(1  .+ ur .^2),

sqrt.(1  .+ ur .^2) .*tauDiffusion,

(etaM .*ur .*( .-6  .+ r .^4  .+ tau .^4)) ./(3. .*sqrt.(1  .+ ur .^2)),

(etaM .*ur .*(3  .- 2 .*r .^4  .+ tau .^4)) ./(3. .*sqrt.(1  .+ ur .^2)),

 .-0.3333333333333333 .*(etaM .*ur .*(3  .+ r .^4  .+ tau .^4)) ./sqrt.(1  .+ ur .^2),

0,

cM .*ur .*sqrt.(1  .+ ur .^2),

(sqrt.(1  .+ ur .^2) .*tauM .*(2  .+ tau .^4)) ./3.,

(sqrt.(1  .+ ur .^2) .*tauM .*( .-1  .+ tau .^4)) ./3.,

 .-0.3333333333333333 .*(sqrt.(1  .+ ur .^2) .*tauM .*( .-1  .+ tau .^4)),

0,

0,

 .-0.3333333333333333 .*(sqrt.(1  .+ ur .^2) .*tauM .*(r .^4  .- tau .^4)),

(sqrt.(1  .+ ur .^2) .*tauM .*(2 .*r .^4  .+ tau .^4)) ./3.,

(sqrt.(1  .+ ur .^2) .*tauM .*(r .^4  .- tau .^4)) ./3.,

0,

cM .*ur .*sqrt.(1  .+ ur .^2),

 .-0.3333333333333333 .*(sqrt.(1  .+ ur .^2) .*tauM .*( .-2  .+ r .^4  .+ tau .^4)),

(sqrt.(1  .+ ur .^2) .*tauM .*( .-1  .+ 2 .*r .^4  .- tau .^4)) ./3.,

(sqrt.(1  .+ ur .^2) .*tauM .*(1  .+ r .^4  .+ tau .^4)) ./3.)
    )
        #########################################################################################################################################################################
        
    Ax=SMatrix{5,5}(
(ur .*dn_dmu,

kappa .*(1  .+ ur .^2),

0,

0,

0,

1,

ur .*tauDiffusion,

(etaM .*( .-6  .+ r .^4  .+ tau .^4)) ./3.,

(etaM .*(3  .- 2 .*r .^4  .+ tau .^4)) ./3.,

 .-0.3333333333333333 .*(etaM .*(3  .+ r .^4  .+ tau .^4)),

0,

cM .*(1  .+ ur .^2),

(ur .*tauM .*(2  .+ tau .^4)) ./3.,

(ur .*tauM .*( .-1  .+ tau .^4)) ./3.,

 .-0.3333333333333333 .*(ur .*tauM .*( .-1  .+ tau .^4)),

0,

0,

(ur .*tauM .*( .-r .^4  .+ tau .^4)) ./3.,

(ur .*tauM .*(2 .*r .^4  .+ tau .^4)) ./3.,

(ur .*tauM .*(r .^4  .- tau .^4)) ./3.,

0,

cM .*(1  .+ ur .^2),

 .-0.3333333333333333 .*(ur .*tauM .*( .-2  .+ r .^4  .+ tau .^4)),

(ur .*tauM .*( .-1  .+ 2 .*r .^4  .- tau .^4)) ./3.,

(ur .*tauM .*(1  .+ r .^4  .+ tau .^4)) ./3.)
    )
        #########################################################################################################################################################################
    
    source=SVector{5}(
(n .*(drur  .+ (dtur .*ur) ./sqrt.(1  .+ ur .^2)  .+ ur ./r  .+ sqrt.(1  .+ ur .^2) ./tau)  .+ X[2] .*(1 ./r  .+ (ur  .+ ur .^3  .+ dtur .*tau) ./((1  .+ ur .^2) .^1.5 .*tau))  .+ (drT .*ur  .+ dtT .*sqrt.(1  .+ ur .^2)) .*dn_dT,

(1  .- (dtur .*ur .*tauDiffusion) ./sqrt.(1  .+ ur .^2)  .+ drur .*( .-1  .+ 1 ./(1  .+ ur .^2)) .*tauDiffusion) .*X[2]  .+ cM .*(drur .*ur  .+ dtur .*sqrt.(1  .+ ur .^2)) .*(X[3]  .+ X[5]),

(3 .*(1  .+ ur .^2) .*(2 .*sqrt.(1  .+ ur .^2) .*X[3]  .+ 2 .*ur .*sqrt.(1  .+ ur .^2) .*tauM .*r .^3 .*X[4]  .- sqrt.(1  .+ ur .^2) .*r .^4 .*X[4]  .+ 2 .*sqrt.(1  .+ ur .^2) .*X[5]  .+ 2 .*ur .*sqrt.(1  .+ ur .^2) .*tauM .*r .^3 .*X[5]  .- sqrt.(1  .+ ur .^2) .*r .^4 .*X[5]  .- 2 .*(1  .+ ur .^2) .*tauM .*(X[3]  .+ X[4]  .- X[5]) .*tau .^3  .+ sqrt.(1  .+ ur .^2) .*X[3] .*tau .^4  .+ sqrt.(1  .+ ur .^2) .*X[4] .*tau .^4  .- sqrt.(1  .+ ur .^2) .*X[5] .*tau .^4)  .+ etaM .*X[2] .*(5 .*dtur .*( .-2  .+ r .^4  .+ tau .^4)  .+ 2 .*dtur .*ur .^2 .*(4  .+ r .^4  .+ tau .^4)  .+ 2 .*drur .*ur .*sqrt.(1  .+ ur .^2) .*(4  .+ r .^4  .+ tau .^4))) ./(9. .*(1  .+ ur .^2) .^1.5),

( .-3 .*(1  .+ ur .^2) .*(sqrt.(1  .+ ur .^2) .*X[3]  .+ 4 .*ur .*sqrt.(1  .+ ur .^2) .*tauM .*r .^3 .*X[4]  .- 2 .*sqrt.(1  .+ ur .^2) .*r .^4 .*X[4]  .+ sqrt.(1  .+ ur .^2) .*X[5]  .+ 4 .*ur .*sqrt.(1  .+ ur .^2) .*tauM .*r .^3 .*X[5]  .- 2 .*sqrt.(1  .+ ur .^2) .*r .^4 .*X[5]  .+ 2 .*(1  .+ ur .^2) .*tauM .*(X[3]  .+ X[4]  .- X[5]) .*tau .^3  .- sqrt.(1  .+ ur .^2) .*X[3] .*tau .^4  .- sqrt.(1  .+ ur .^2) .*X[4] .*tau .^4  .+ sqrt.(1  .+ ur .^2) .*X[5] .*tau .^4)  .+ etaM .*X[2] .*(2 .*dtur .*ur .^2 .*( .-2  .- 2 .*r .^4  .+ tau .^4)  .+ 2 .*drur .*ur .*sqrt.(1  .+ ur .^2) .*( .-2  .- 2 .*r .^4  .+ tau .^4)  .+ 5 .*dtur .*(1  .- 2 .*r .^4  .+ tau .^4))) ./(9. .*(1  .+ ur .^2) .^1.5),

(3 .*(1  .+ ur .^2) .*(sqrt.(1  .+ ur .^2) .*X[3]  .- 2 .*ur .*sqrt.(1  .+ ur .^2) .*tauM .*r .^3 .*X[4]  .+ sqrt.(1  .+ ur .^2) .*r .^4 .*X[4]  .+ sqrt.(1  .+ ur .^2) .*X[5]  .- 2 .*ur .*sqrt.(1  .+ ur .^2) .*tauM .*r .^3 .*X[5]  .+ sqrt.(1  .+ ur .^2) .*r .^4 .*X[5]  .+ 2 .*(1  .+ ur .^2) .*tauM .*(X[3]  .+ X[4]  .- X[5]) .*tau .^3  .- sqrt.(1  .+ ur .^2) .*X[3] .*tau .^4  .- sqrt.(1  .+ ur .^2) .*X[4] .*tau .^4  .+ sqrt.(1  .+ ur .^2) .*X[5] .*tau .^4)  .- etaM .*X[2] .*(2 .*dtur .*ur .^2 .*( .-2  .+ r .^4  .+ tau .^4)  .+ 2 .*drur .*ur .*sqrt.(1  .+ ur .^2) .*( .-2  .+ r .^4  .+ tau .^4)  .+ 5 .*dtur .*(1  .+ r .^4  .+ tau .^4))) ./(9. .*(1  .+ ur .^2) .^1.5))
    )
    return (At,Ax, source) 
end

    