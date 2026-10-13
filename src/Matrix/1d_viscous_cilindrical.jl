# 1+1D cylindrical viscous hydro: Israel-Stewart + the DNMR term delta_pipi pi theta (switchable).
# Fields phi = (T, u^r, pi^phi_phi, pi^eta_eta, Pi)  (the old line here said "T, ux, pi^yy, pi^zz, piB").
#   Nabla_m T^{mn} = 0;  tauS Delta D pi + pi + 2 eta sigma + dpp pi theta = 0;  tauB D Pi + Pi + zeta theta = 0.
#
# delta_pipi (2026-10-13, RS: "4/3 tau_pi everywhere"; MATRICES_PROGRAMME.md P2b). `dpp` is the last
# argument of `one_d_viscous_matrix` and defaults to `delta_pipi(tauS)` = FLUIDUM_DELTA_PIPI * (4/3) tauS:
#   FLUIDUM_DELTA_PIPI unset or "1" (DEFAULT) -> dpp = 4/3 tau_pi, the medium of the 7-field :MIS kernel
#                                               (1d_viscous_HQ_cilindrical_fugacity.jl), exactly;
#   FLUIDUM_DELTA_PIPI = "0"                  -> dpp = 0, plain IS: bit-for-bit what this kernel gave
#                                               before 2026-10-13 (the dpp = 0 branch returns early).
# Until 2026-10-13 this kernel was plain IS (delta_pipi = 0) while the :MIS kernel carried 4/3 tau_pi.
# The term is stored as its own coefficient matrices (dAt, dAx, dsource), so At, Ax, source below are
# literally the pre-10-13 code.
#
# GATED against the covariant derivation, read out of THIS file: gate M-V5 of
# Julia/tools/derive_1d_matrices.wls (raw A_t, A_r, src, symbolically, for symbolic dpp; at
# dpp = 4/3 tauS the medium block of the 7-field :MIS kernel) and V-G4 of derive_2p1d_viscous.wls;
# J-V5/J-V5d of Projects/FluidumValidation/gate_fluidum_matrices.jl (programme.jl check, both switch
# values). Also D1 of bench/bench_1p1d_viscous.jl and B-G1 of
# Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl against `one_d_viscous_matrix_derived`.



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



# ⛔ RETRACTED IN PLACE 2026-10-13: `matrix1d_visc` (the non-mutating twin) cannot run and never
# could on any shipped EoS -- `pressure(ϕ[1], Val(1), params.eos)` has no method (MethodError:
# thermodynamic(::Float64, ::Val{1}, ::FluiduMEoS)); the intended call is `pressure(ϕ[1], params.eos)`.
# Nothing in the package, Julia/Projects or Tex calls it. Kept, unrepaired, because this repo keeps
# its wrong turns; use `matrix1d_visc!` above.
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
    # theta = (u^r/u^tau) d_tau u^r + d_r u^r + u^tau/tau + u^r/r  (gate M-V5, symbolic dpp)
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
twin and the 2-D medium (`two_d_viscous_matrix_derived`): 4/3 tau_pi by default (RS, 2026-10-13),
0 when `FLUIDUM_DELTA_PIPI=0` (plain IS, pre-2026-10-13 results bit-for-bit). `DELTA_PIPI` is read in
`__init__`; assign `Fluidum.DELTA_PIPI[] = 0.0` in-process to switch. The 7-field `:MIS` kernel
(`one_d_viscous_matrix_fugacity`) carries 4/3 tau_pi unconditionally and ignores the switch; the
retracted legacy `two_d_viscous_matrix` carries none.
"""
@inline delta_pipi(tauS) = DELTA_PIPI[] * (4/3) * tauS


