# 1+1D cylindrical viscous hydro: Israel-Stewart + the DNMR term delta_pipi pi theta (switchable).
# Fields phi = (T, u^r, pi^phi_phi, pi^eta_eta, Pi).
#   Nabla_m T^{mn} = 0;  tauS Delta D pi + pi + 2 eta sigma + dpp pi theta = 0;  tauB D Pi + Pi + zeta theta = 0.
#
# `dpp` is the last argument of `one_d_viscous_matrix` and defaults to `delta_pipi(tauS)`
# = FLUIDUM_DELTA_PIPI * (4/3) tauS: 4/3 tau_pi by default (the medium of the 7-field :MIS kernel,
# 1d_viscous_HQ_cilindrical_fugacity.jl), 0 for plain IS with FLUIDUM_DELTA_PIPI=0. The term is stored
# as separate coefficient matrices (dAt, dAx, dsource). The matrices agree with the covariant
# derivation and with `one_d_viscous_matrix_derived` (deriv/derive_2p1d_viscous.wls).



function matrix1d_visc!(A_i,Source,ϕ,t,X,params)
 
    p = pressure(ϕ[1],params.eos) #entropy
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
    
 
    (At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],p,dpt,dptt,zeta,etaVisc,tauS,tauB)
    

        
        Ainv= inv(At)

    
        A_mul_B!(A_i[1], Ainv,Ax)
        
        
        jgemvavx!(Source, Ainv,source)
    
        
    
    end 



# Not usable: `pressure(ϕ[1], Val(1), params.eos)` has no method, so `matrix1d_visc` throws a
# MethodError on every EoS. Use `matrix1d_visc!` above.
function matrix1d_visc(ϕ,t,X,params)
 

    p = pressure(ϕ[1],Val(1),params.eos) #entropy
    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
            
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
        
     
        (At,Ax, source)=one_d_viscous_matrix(ϕ,t,X[1],p,dpt,dptt,zeta,etaVisc,tauS,tauB)
        
    
            
        Ainv= inv(At)
    
        
          
        return (Ainv*source,Ainv*Ax)
            
        
end 



@inbounds function one_d_viscous_matrix(u,tau,R,p,dtp,dtdtp,zeta,visc,tauS,tauB,dpp=delta_pipi(tauS))
    
At=SMatrix{5,5}(

    dtdtp*u[1] + (dtp + dtdtp*u[1])*^(u[2],2)

    ,(dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,0

    ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])

    ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))

    ,(-2*u[2]*visc)/(3*^(R,2)*sqrt(1 + ^(u[2],2)))

    ,(-2*u[2]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2)))

    ,(u[2]*zeta)/sqrt(1 + ^(u[2],2))

    ,-^(u[2],2)

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,(tauS*sqrt(1 + ^(u[2],2)))/^(R,2)

    ,0

    ,0

    ,-^(u[2],2)

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,0

    ,(tauS*sqrt(1 + ^(u[2],2)))/^(tau,2)

    ,0

    ,^(u[2],2)

    ,u[2]*sqrt(1 + ^(u[2],2))

    ,0

    ,0

    ,tauB*sqrt(1 + ^(u[2],2))

)
    #########################################################################################################################################################################
    
Ax=SMatrix{5,5}(

    (dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2))

    ,dtp + (dtp + dtdtp*u[1])*^(u[2],2)

    ,0

    ,0

    ,0

    ,((1 + 2*^(u[2],2))*(dtp*u[1] + u[5] - u[3] - u[4]))/sqrt(1 + ^(u[2],2))

    ,2*u[2]*(dtp*u[1] + u[5] - u[3] - u[4])

    ,(-2*visc)/(3*^(R,2))

    ,(-2*visc)/(3*^(tau,2))

    ,zeta

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,-1 - ^(u[2],2)

    ,(tauS*u[2])/^(R,2)

    ,0

    ,0

    ,-(u[2]*sqrt(1 + ^(u[2],2)))

    ,-1 - ^(u[2],2)

    ,0

    ,(tauS*u[2])/^(tau,2)

    ,0

    ,u[2]*sqrt(1 + ^(u[2],2))

    ,1 + ^(u[2],2)

    ,0

    ,0

    ,tauB*u[2]

)
    #########################################################################################################################################################################

source=SVector{5}(

    (dtp*u[1]*(R + R*^(u[2],2) + tau*u[2]*sqrt(1 + ^(u[2],2))) + tau*u[2]*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + R*((1 + ^(u[2],2))*u[5] + u[4] - ^(u[2],2)*(u[3] + u[4])))/(R*tau)

    ,(dtp*u[1]*u[2]*(tau*u[2] + R*sqrt(1 + ^(u[2],2))) + R*u[2]*sqrt(1 + ^(u[2],2))*(u[5] - u[3] - u[4]) + tau*(-2*u[3] + ^(u[2],2)*(u[5] - u[3] - u[4]) - u[4]))/(R*tau)

    ,(4*tau*u[2]*visc - 2*R*sqrt(1 + ^(u[2],2))*visc + 3*R*tau*u[3])/(3*^(R,3)*tau)

    ,(-2*tau*u[2]*visc + 4*R*sqrt(1 + ^(u[2],2))*visc + 3*R*tau*u[4])/(3*R*^(tau,3))

    ,(u[2]*zeta)/R + (sqrt(1 + ^(u[2],2))*zeta)/tau + u[5]

)

iszero(dpp) && return (At,Ax, source)
    # delta_pipi: the coefficient of dpp of  dpp pi theta  in the two shear rows,
    # theta = (u^r/u^tau) d_tau u^r + d_r u^r + u^tau/tau + u^r/r
dAt=SMatrix{5,5}(
    0,0,0,0,0
    ,0,0,(u[3]*u[2])/(^(R,2)*sqrt(1 + ^(u[2],2))),(u[4]*u[2])/(^(tau,2)*sqrt(1 + ^(u[2],2))),0
    ,0,0,0,0,0
    ,0,0,0,0,0
    ,0,0,0,0,0
)
dAx=SMatrix{5,5}(
    0,0,0,0,0
    ,0,0,u[3]/^(R,2),u[4]/^(tau,2),0
    ,0,0,0,0,0
    ,0,0,0,0,0
    ,0,0,0,0,0
)
dsource=SVector{5}(
    0,0
    ,(u[3]*(sqrt(1 + ^(u[2],2))/tau + u[2]/R))/^(R,2)
    ,(u[4]*(sqrt(1 + ^(u[2],2))/tau + u[2]/R))/^(tau,2)
    ,0
)
return (At + dpp*dAt, Ax + dpp*dAx, source + dpp*dsource)
end

"""
    delta_pipi(tauS) = DELTA_PIPI[] * (4/3) * tauS

The DNMR coefficient delta_pipi used by the 5-field 1-D kernel (`one_d_viscous_matrix`), its generated
counterpart and the 2-D medium (`two_d_viscous_matrix_derived`): 4/3 tau_pi by default, 0 (plain IS)
when `FLUIDUM_DELTA_PIPI=0`. `DELTA_PIPI` is read in `__init__`; assign `Fluidum.DELTA_PIPI[] = 0.0`
in-process to switch. The 7-field `:MIS` kernel (`one_d_viscous_matrix_fugacity`) always uses
4/3 tau_pi; the legacy `two_d_viscous_matrix` has none.
"""
@inline delta_pipi(tauS) = DELTA_PIPI[] * (4/3) * tauS


