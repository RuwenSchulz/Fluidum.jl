# the convention here are T, ux, uy, \[Pi]yy, \[Pi]zz, \[Pi]xy, \[Pi]B
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
# ⛔ 2026-09-04 — DO NOT USE.  RETRACTED IN PLACE, REPLACED BY Matrix/HQ_2p1d_BG.jl.
#
# `matrix2d_visc_HQ!` and its `two_d_viscous_HQ_matrix` are the ORIGINAL 10-field 2+1D
# hydro+charm system, and they are unusable — the note recording that lived only in
# HQ_2p1d_BG.jl until today, while this file opened with a one-line convention comment as if
# it were live.  What is wrong, measured (Projects/FluidumValidation):
#
#   * Characteristic speeds: max|Im λ| = 3.02, max|Re λ| = 36.5 c over the physical state
#     lattice of bench_fluidum_full.jl (98–130 c with complex pairs on the harder probe of
#     bench_fluidum_2p1d.jl).  This is the SOLE entry behind that audit's standing M3
#     SUPERLUMINAL failure and one of four behind M2.  The replacement `matrix2d_visc_HQ_BG!`
#     measures max|Im λ| = 3.1e-16, max|Re λ| = 0.962 on the same lattice.
#   * The wrapper reads `ϕ[6]` (= π^{xy} in its own documented layout) as the fugacity.
#   * The wrapper passes `dpt` (entropy) into the matrix's `p` (pressure) slot — inert only
#     because the body never uses `p`.
#   * Its hydro block predates the 2026-09-03 shear-row repair of Matrix/2d_viscous.jl and is
#     NOT routed through `two_d_viscous_matrix_active`, so it carries the ten wrong entries
#     unconditionally.
#
# Kept unreferenced rather than deleted — this repo preserves wrong turns, dated.  Use
# `matrix2d_visc_HQ_BG!` (HQ_2p1d_BG.jl), whose hydro block is `matrix2d_visc!`'s and whose
# charm rows are the derived `hq2d_live_rows`, gated L1–L4 in bench_fluidum_2p1d.jl.
# ─────────────────────────────────────────────────────────────────────────────────────────────

@inbounds @fastmath function matrix2d_visc_HQ!(A_i,Source,ϕ,t,X,params)

    dpt = pressure_derivative(ϕ[1],Val(1),params.eos) #entropy
    dptt = pressure_derivative(ϕ[1],Val(2),params.eos)
        
    etaVisc=viscosity(ϕ[1],dpt,params.shear)
    tauS=τ_shear(ϕ[1],dpt,params.shear)
    tauB=τ_bulk(ϕ[1],dpt,dptt,params.bulk)
    zeta=bulk_viscosity(ϕ[1],dpt,params.bulk)
   
    thermo = thermodynamic(ϕ[1],ϕ[6],params.eos.hadron_list)
    n=thermo.pressure
    dtn, dmn = thermo.pressure_derivative
    dmn+=0.0001
    Ds = diffusion(ϕ[1],n,params.diffusion)
    tauDiff=τ_diffusion(ϕ[1],params.diffusion)

    #κ = diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #diffusion coefficient for hadrons
    #tauDiff=τ_diffusion_hadron(ϕ[1],ϕ[6],params.eos,params.diffusion) #tau diffusion for hadrons
    
    dmp = 0 #for now we don t have chemical potential in the eos
    dtdmp = 0 
    dmdmp = 0
 
    (At,Ax, Ay, source)=two_d_viscous_HQ_matrix(ϕ,t,dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    

        
        Ainv= inv(At)

    
        A_mul_B!(A_i[1], Ainv,Ax)
        A_mul_B!(A_i[2], Ainv,Ay)
        
        @inbounds for i in 1:7
          @inbounds for j in 8:10
            A_i[1][i,j] = 0.
            A_i[2][i,j] = 0.
            end
        end

        jgemvavx!(Source, Ainv,source)
    
        
    
end 



@inbounds @fastmath function two_d_viscous_HQ_matrix(u,tau,p,dtp,dtdtp,zeta,visc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)

    #u is a vector of the form (T, ux, uy, \[Pi]yy, \[Pi]zz, \[Pi]xy, \[Pi]B, α, nux, nuy)
    #tau is the time step
    #p is the pressure
    #dtp is the first derivative of the pressure with respect to T
    #dtdtp is the second derivative of the pressure with respect to T
    #zeta is the bulk viscosity
    #visc is the shear viscosity
    #tauS is the shear relaxation time
    #tauB is the bulk relaxation time
    #n is the number density
    #dtn is the first derivative of n with respect to T
    #dmn is the second derivative of n with respect to T
    #tauDiff is the diffusion relaxation time
    #Ds is the diffusion coefficient
   
At=SMatrix{10,10}(

   dtp*(^(getindex(u,2),2) + ^(getindex(u,3),2)) + dtdtp*getindex(u,1)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(dtp + dtdtp*getindex(u,1))*getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(dtp + dtdtp*getindex(u,1))*getindex(u,3)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,0

    ,dtn*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,2*(dtp*getindex(u,1)*getindex(u,2) + (-(getindex(u,2)*(getindex(u,4) + getindex(u,5))) + getindex(u,3)*getindex(u,6))/(1 + ^(getindex(u,3),2)) + getindex(u,2)*getindex(u,7))

    ,(dtp*getindex(u,1)*(1 + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + 2*^(getindex(u,2),2) + ^(getindex(u,3),2)) - getindex(u,4) - getindex(u,5) + 2*^(getindex(u,2),3)*getindex(u,3)*getindex(u,6) + 3*getindex(u,2)*getindex(u,3)*(1 + ^(getindex(u,3),2))*getindex(u,6) + getindex(u,7) - 2*^(getindex(u,2),4)*(getindex(u,4) + getindex(u,5) - (1 + ^(getindex(u,3),2))*getindex(u,7)) + 3*^(getindex(u,2),2)*(1 + ^(getindex(u,3),2))*(-getindex(u,4) - getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7)) + ^(getindex(u,3),2)*(-getindex(u,4) - (2 + ^(getindex(u,3),2))*getindex(u,5) + (3 + 3*^(getindex(u,3),2) + ^(getindex(u,3),4))*getindex(u,7)))/((1 + ^(getindex(u,3),2))*^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5))

    ,(getindex(u,6) + getindex(u,3)*(dtp*getindex(u,1)*getindex(u,2)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + getindex(u,3)*getindex(u,6) + getindex(u,2)*(-getindex(u,4) + (1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*getindex(u,7))))/^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5)

    ,(-2*(getindex(u,2)*(visc + ^(getindex(u,3),2)*(visc - 3*tauS*getindex(u,4))) + 3*tauS*getindex(u,3)*(1 + ^(getindex(u,3),2))*getindex(u,6)))/(3. *sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(-2*visc*getindex(u,2))/(3. *^(tau,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(getindex(u,3)*((3 + 6*tauS)*getindex(u,4) + 6*tauS*(1 + ^(getindex(u,3),2))*getindex(u,5) + ^(getindex(u,2),2)*(-4*visc + 6*tauS*(2*getindex(u,4) + getindex(u,5)))) - 3*getindex(u,2)*(-1 + 2*tauS + 4*tauS*^(getindex(u,3),2))*getindex(u,6))/(6. *sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(zeta*getindex(u,2))/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(n*^(getindex(u,2),3) + (1 + ^(getindex(u,3),2))*getindex(u,9) + getindex(u,2)*(n + n*^(getindex(u,3),2) - getindex(u,3)*getindex(u,10)))/^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5)

    ,(tauDiff*getindex(u,2)*(-((1 + ^(getindex(u,3),2))*getindex(u,9)) + getindex(u,2)*getindex(u,3)*getindex(u,10)))/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,-((tauDiff*getindex(u,3)*((1 + ^(getindex(u,3),2))*getindex(u,9) - getindex(u,2)*getindex(u,3)*getindex(u,10)))/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(2*(dtp*getindex(u,1)*getindex(u,3)*^(1 + ^(getindex(u,3),2),2) + getindex(u,2)*getindex(u,6) - getindex(u,2)*^(getindex(u,3),2)*getindex(u,6) + 2*^(getindex(u,3),3)*getindex(u,7) + ^(getindex(u,3),5)*getindex(u,7) + getindex(u,3)*((1 + ^(getindex(u,2),2))*getindex(u,4) + ^(getindex(u,2),2)*getindex(u,5) + getindex(u,7))))/^(1 + ^(getindex(u,3),2),2)

    ,(dtp*getindex(u,1)*getindex(u,2)*getindex(u,3)*^(1 + ^(getindex(u,3),2),2)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + 2*^(getindex(u,2),5)*getindex(u,3)*(getindex(u,4) + getindex(u,5)) - 2*^(getindex(u,2),4)*(-1 + ^(getindex(u,3),2))*getindex(u,6) + ^(1 + ^(getindex(u,3),2),2)*getindex(u,6) - 3*^(getindex(u,2),2)*(-1 + ^(getindex(u,3),4))*getindex(u,6) + getindex(u,2)*getindex(u,3)*(1 + ^(getindex(u,3),2))*(3*getindex(u,4) + (1 + ^(getindex(u,3),2))*(getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7))) + ^(getindex(u,2),3)*getindex(u,3)*((5 + 3*^(getindex(u,3),2))*getindex(u,4) + (1 + ^(getindex(u,3),2))*(3*getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7))))/(^(1 + ^(getindex(u,3),2),2)*^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5))

    ,(dtp*getindex(u,1)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + 2*^(getindex(u,3),2)) + (1 + ^(getindex(u,2),2))*getindex(u,4) - getindex(u,2)*getindex(u,3)*getindex(u,6) + (1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + 2*^(getindex(u,3),2))*getindex(u,7))/^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5)

    ,(getindex(u,3)*(4*visc*(1 + ^(getindex(u,3),2)) + 3*(1 - 2*tauS)*getindex(u,4) + 6*^(getindex(u,2),2)*(visc - tauS*getindex(u,4))) + 3*getindex(u,2)*(1 + 2*tauS*^(getindex(u,3),2))*getindex(u,6))/(3. *sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(-2*visc*getindex(u,3))/(3. *^(tau,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(getindex(u,2)*(-3*(-4*visc + getindex(u,4) + 2*tauS*getindex(u,4) + getindex(u,5)) + ^(getindex(u,3),4)*(8*visc - 6*tauS*getindex(u,5)) + ^(getindex(u,3),2)*(20*visc - 3*getindex(u,5) - 6*tauS*(2*getindex(u,4) + getindex(u,5))) - 3*^(getindex(u,2),2)*(-4*visc + getindex(u,4) + 2*tauS*getindex(u,4) + getindex(u,5) + ^(getindex(u,3),2)*(-4*visc + 4*tauS*getindex(u,4) + 2*tauS*getindex(u,5)))) + 3*getindex(u,3)*(-((-1 + 2*tauS)*(1 + ^(getindex(u,3),2))) + ^(getindex(u,2),2)*(2 + 4*tauS*^(getindex(u,3),2)))*getindex(u,6))/(6. *(1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(zeta*getindex(u,3))/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(getindex(u,3)*(n*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) - getindex(u,2)*getindex(u,9)) + (1 + ^(getindex(u,2),2))*getindex(u,10))/^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5)

    ,(tauDiff*^(getindex(u,2),2)*getindex(u,3)*getindex(u,9) - tauDiff*(getindex(u,2) + ^(getindex(u,2),3))*getindex(u,10))/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(tauDiff*getindex(u,3)*(getindex(u,2)*getindex(u,3)*getindex(u,9) - (1 + ^(getindex(u,2),2))*getindex(u,10)))/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,1 - (1 + ^(getindex(u,2),2))/(1 + ^(getindex(u,3),2))

    ,-((getindex(u,2)*(1 + ^(getindex(u,2),2)))/((1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))))

    ,getindex(u,3)/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,tauS*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,-(^(getindex(u,2),2)/(1 + ^(getindex(u,3),2)))

    ,-((getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))/(1 + ^(getindex(u,3),2)))

    ,0

    ,0

    ,(tauS*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))/^(tau,2)

    ,0

    ,0

    ,0

    ,0

    ,0

    ,(2*getindex(u,2)*getindex(u,3))/(1 + ^(getindex(u,3),2))

    ,(getindex(u,3)*(1 + 2*^(getindex(u,2),2) + ^(getindex(u,3),2)))/((1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,getindex(u,2)/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,tauS*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,0

    ,^(getindex(u,2),2) + ^(getindex(u,3),2)

    ,getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,getindex(u,3)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,tauB*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,dmn*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,Ds*n*getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,Ds*n*getindex(u,3)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,getindex(u,2)/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,tauDiff*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,getindex(u,3)/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,tauDiff*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))
)
    #########################################################################################################################################################################
    
Ax=SMatrix{10,10}(
    (dtp + dtdtp*getindex(u,1))*getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,dtp + (dtp + dtdtp*getindex(u,1))*^(getindex(u,2),2)

    ,(dtp + dtdtp*getindex(u,1))*getindex(u,2)*getindex(u,3)

    ,0

    ,0

    ,0

    ,0

    ,dtn*getindex(u,2)

    ,0

    ,0

    ,(dtp*getindex(u,1)*(1 + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + 2*^(getindex(u,2),2) + ^(getindex(u,3),2)) - getindex(u,4) - getindex(u,5) + 2*^(getindex(u,2),3)*getindex(u,3)*getindex(u,6) + 3*getindex(u,2)*getindex(u,3)*(1 + ^(getindex(u,3),2))*getindex(u,6) + getindex(u,7) - 2*^(getindex(u,2),4)*(getindex(u,4) + getindex(u,5) - (1 + ^(getindex(u,3),2))*getindex(u,7)) + 3*^(getindex(u,2),2)*(1 + ^(getindex(u,3),2))*(-getindex(u,4) - getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7)) + ^(getindex(u,3),2)*(-getindex(u,4) - (2 + ^(getindex(u,3),2))*getindex(u,5) + (3 + 3*^(getindex(u,3),2) + ^(getindex(u,3),4))*getindex(u,7)))/((1 + ^(getindex(u,3),2))*^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5))

    ,2*(dtp*getindex(u,1)*getindex(u,2) + (-(getindex(u,2)*(getindex(u,4) + getindex(u,5))) + getindex(u,3)*getindex(u,6))/(1 + ^(getindex(u,3),2)) + getindex(u,2)*getindex(u,7))

    ,getindex(u,3)*(dtp*getindex(u,1) + getindex(u,7))

    ,(-2*(visc*^(1 + ^(getindex(u,3),2),2) + ^(getindex(u,2),2)*(visc + ^(getindex(u,3),2)*(visc - 3*tauS*getindex(u,4))) + 3*tauS*getindex(u,2)*getindex(u,3)*(1 + ^(getindex(u,3),2))*getindex(u,6)))/(3. *(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(-2*visc)/(3. *^(tau,2))

    ,(getindex(u,2)*(getindex(u,3)*(-4*visc - 4*visc*^(getindex(u,3),2) + (3 + 6*tauS)*getindex(u,4) - 4*^(getindex(u,2),2)*(visc - 3*tauS*getindex(u,4)) + 6*tauS*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*getindex(u,5)) - 3*getindex(u,2)*(-1 + 2*tauS + 4*tauS*^(getindex(u,3),2))*getindex(u,6)))/(6. *(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,zeta

    ,n

    ,(tauDiff*^(getindex(u,2),2)*(-((1 + ^(getindex(u,3),2))*getindex(u,9)) + getindex(u,2)*getindex(u,3)*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(tauDiff*getindex(u,2)*getindex(u,3)*(-((1 + ^(getindex(u,3),2))*getindex(u,9)) + getindex(u,2)*getindex(u,3)*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(dtp*getindex(u,1)*getindex(u,2)*getindex(u,3)*^(1 + ^(getindex(u,3),2),2)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + 2*^(getindex(u,2),5)*getindex(u,3)*(getindex(u,4) + getindex(u,5)) - 2*^(getindex(u,2),4)*(-1 + ^(getindex(u,3),2))*getindex(u,6) + ^(1 + ^(getindex(u,3),2),2)*getindex(u,6) - 3*^(getindex(u,2),2)*(-1 + ^(getindex(u,3),4))*getindex(u,6) + getindex(u,2)*getindex(u,3)*(1 + ^(getindex(u,3),2))*(3*getindex(u,4) + (1 + ^(getindex(u,3),2))*(getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7))) + ^(getindex(u,2),3)*getindex(u,3)*((5 + 3*^(getindex(u,3),2))*getindex(u,4) + (1 + ^(getindex(u,3),2))*(3*getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7))))/(^(1 + ^(getindex(u,3),2),2)*^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5))

    ,(2*getindex(u,3)*((1 + ^(getindex(u,2),2))*getindex(u,4) + ^(getindex(u,2),2)*getindex(u,5)) - 2*getindex(u,2)*(-1 + ^(getindex(u,3),2))*getindex(u,6))/^(1 + ^(getindex(u,3),2),2)

    ,getindex(u,2)*(dtp*getindex(u,1) + getindex(u,7))

    ,getindex(u,6) + (2*getindex(u,2)*getindex(u,3)*(visc + visc*^(getindex(u,3),2) - tauS*getindex(u,4) + ^(getindex(u,2),2)*(visc - tauS*getindex(u,4)) + tauS*getindex(u,2)*getindex(u,3)*getindex(u,6)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,0

    ,((1 + ^(getindex(u,3),2))*(4*visc*(1 + ^(getindex(u,3),2)) - 2*getindex(u,4) - (1 + ^(getindex(u,3),2))*getindex(u,5)) + ^(getindex(u,2),2)*(8*visc - (3 + 2*tauS)*getindex(u,4) - 2*getindex(u,5) + ^(getindex(u,3),4)*(4*visc - 2*tauS*getindex(u,5)) - 2*^(getindex(u,3),2)*(-6*visc + getindex(u,4) + 2*tauS*getindex(u,4) + (1 + tauS)*getindex(u,5))) + ^(getindex(u,2),4)*(4*visc - (1 + 2*tauS)*getindex(u,4) - getindex(u,5) + ^(getindex(u,3),2)*(4*visc - 2*tauS*(2*getindex(u,4) + getindex(u,5)))) - (-3 + 2*tauS)*getindex(u,2)*getindex(u,3)*(1 + ^(getindex(u,3),2))*getindex(u,6) + 2*^(getindex(u,2),3)*getindex(u,3)*(1 + 2*tauS*^(getindex(u,3),2))*getindex(u,6))/(2. *(1 + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,0

    ,0

    ,(tauDiff*^(getindex(u,2),2)*(getindex(u,2)*getindex(u,3)*getindex(u,9) - (1 + ^(getindex(u,2),2))*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,(tauDiff*getindex(u,3)*(^(getindex(u,2),2)*getindex(u,3)*getindex(u,9) - (getindex(u,2) + ^(getindex(u,2),3))*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,-((getindex(u,2)*(1 + ^(getindex(u,2),2)))/((1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))))

    ,-((1 + ^(getindex(u,2),2))/(1 + ^(getindex(u,3),2)))

    ,0

    ,tauS*getindex(u,2)

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,-((getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))/(1 + ^(getindex(u,3),2)))

    ,-1 - ^(getindex(u,2),2)/(1 + ^(getindex(u,3),2))

    ,0

    ,0

    ,(tauS*getindex(u,2))/^(tau,2)

    ,0

    ,0

    ,0

    ,0

    ,0

    ,(getindex(u,3)*(1 + 2*^(getindex(u,2),2) + ^(getindex(u,3),2)))/((1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

    ,(2*getindex(u,2)*getindex(u,3))/(1 + ^(getindex(u,3),2))

    ,1

    ,0

    ,0

    ,tauS*getindex(u,2)

    ,0

    ,0

    ,0

    ,0

    ,getindex(u,2)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

    ,1 + ^(getindex(u,2),2)

    ,getindex(u,2)*getindex(u,3)

    ,0

    ,0

    ,0

    ,tauB*getindex(u,2)

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,dmn*getindex(u,2)

    ,Ds*n*(1 + ^(getindex(u,2),2))

    ,Ds*n*getindex(u,2)*getindex(u,3)

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,1

    ,tauDiff*getindex(u,2)

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,tauDiff*getindex(u,2) 
   
)
    #########################################################################################################################################################################
Ay=SMatrix{10,10}(

(dtp + dtdtp*getindex(u,1))*getindex(u,3)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,(dtp + dtdtp*getindex(u,1))*getindex(u,2)*getindex(u,3)

,dtp + (dtp + dtdtp*getindex(u,1))*^(getindex(u,3),2)

,0

,0

,0

,0

,dtn*getindex(u,3)

,0

,0

,(getindex(u,6) + getindex(u,3)*(dtp*getindex(u,1)*getindex(u,2)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + getindex(u,3)*getindex(u,6) + getindex(u,2)*(-getindex(u,4) + (1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*getindex(u,7))))/^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5)

,getindex(u,3)*(dtp*getindex(u,1) + getindex(u,7))

,0

,-(((1 + 2*tauS*^(getindex(u,3),2))*(-(getindex(u,2)*getindex(u,3)*getindex(u,4)) + (1 + ^(getindex(u,3),2))*getindex(u,6)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,0

,((2 + (1 + 2*tauS)*^(getindex(u,3),2) + ^(getindex(u,2),2)*(2 + 4*tauS*^(getindex(u,3),2)))*getindex(u,4) + (1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + 2*tauS*^(getindex(u,3),2))*getindex(u,5) - getindex(u,2)*getindex(u,3)*(1 + 2*tauS + 4*tauS*^(getindex(u,3),2))*getindex(u,6))/(2. *(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,0

,0

,(tauDiff*getindex(u,2)*getindex(u,3)*(-((1 + ^(getindex(u,3),2))*getindex(u,9)) + getindex(u,2)*getindex(u,3)*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,(tauDiff*^(getindex(u,3),2)*(-((1 + ^(getindex(u,3),2))*getindex(u,9)) + getindex(u,2)*getindex(u,3)*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,(dtp*getindex(u,1)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + 2*^(getindex(u,3),2)) + (1 + ^(getindex(u,2),2))*getindex(u,4) - getindex(u,2)*getindex(u,3)*getindex(u,6) + (1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + 2*^(getindex(u,3),2))*getindex(u,7))/^(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2),1.5)

,getindex(u,2)*(dtp*getindex(u,1) + getindex(u,7))

,2*getindex(u,3)*(dtp*getindex(u,1) + getindex(u,7))

,(4*visc*^(1 + ^(getindex(u,3),2),2) + 3*(1 - 2*tauS)*^(getindex(u,3),2)*getindex(u,4) + ^(getindex(u,2),2)*(4*visc + ^(getindex(u,3),2)*(4*visc - 6*tauS*getindex(u,4))) + 3*getindex(u,2)*getindex(u,3)*(1 + 2*tauS*^(getindex(u,3),2))*getindex(u,6))/(3. *(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,(-2*visc)/(3. *^(tau,2))

,(getindex(u,3)*(getindex(u,2)*(8*visc - 3*(1 + 2*tauS)*getindex(u,4) - 3*getindex(u,5) + ^(getindex(u,3),2)*(8*visc*(2 + ^(getindex(u,3),2)) - 12*tauS*getindex(u,4) - 3*(1 + 2*tauS + 2*tauS*^(getindex(u,3),2))*getindex(u,5)) + ^(getindex(u,2),2)*(8*visc - 3*(1 + 2*tauS)*getindex(u,4) - 3*getindex(u,5) + 2*^(getindex(u,3),2)*(4*visc - 3*tauS*(2*getindex(u,4) + getindex(u,5))))) + 3*getindex(u,3)*(-((-1 + 2*tauS)*(1 + ^(getindex(u,3),2))) + ^(getindex(u,2),2)*(2 + 4*tauS*^(getindex(u,3),2)))*getindex(u,6)))/(6. *(1 + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,zeta

,n

,(tauDiff*getindex(u,3)*(^(getindex(u,2),2)*getindex(u,3)*getindex(u,9) - (getindex(u,2) + ^(getindex(u,2),3))*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,(tauDiff*^(getindex(u,3),2)*(getindex(u,2)*getindex(u,3)*getindex(u,9) - (1 + ^(getindex(u,2),2))*getindex(u,10)))/(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,getindex(u,3)/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,0

,1

,tauS*getindex(u,3)

,0

,0

,0

,0

,0

,0

,0

,0

,0

,0

,(tauS*getindex(u,3))/^(tau,2)

,0

,0

,0

,0

,0

,getindex(u,2)/sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,1

,0

,0

,0

,tauS*getindex(u,3)

,0

,0

,0

,0

,getindex(u,3)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))

,getindex(u,2)*getindex(u,3)

,1 + ^(getindex(u,3),2)

,0

,0

,0

,tauB*getindex(u,3)

,0

,0

,0

,0

,0

,0

,0

,0

,0

,0

,dmn*getindex(u,3)

,Ds*n*getindex(u,2)*getindex(u,3)

,Ds*n*(1 + ^(getindex(u,3),2))

,0

,0

,0

,0

,0

,0

,0

,0

,tauDiff*getindex(u,3)

,0

,0

,0

,0

,0

,0

,0

,0

,1

,0

,tauDiff*getindex(u,3)
   )
    #########################################################################################################################################################


    #########################################################################################################################################################

source=SVector{10}(
(dtp*getindex(u,1)*(1 + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + getindex(u,5) + 2*getindex(u,2)*getindex(u,3)*getindex(u,6) + getindex(u,7) + ^(getindex(u,2),2)*(-getindex(u,4) - getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7)) + ^(getindex(u,3),2)*(getindex(u,4) + getindex(u,5) + (2 + ^(getindex(u,3),2))*getindex(u,7)))/(tau*(1 + ^(getindex(u,3),2)))

,(dtp*getindex(u,1)*getindex(u,2)*(1 + ^(getindex(u,3),2))*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + 2*^(getindex(u,2),2)*getindex(u,3)*getindex(u,6) + getindex(u,3)*(1 + ^(getindex(u,3),2))*getindex(u,6) + ^(getindex(u,2),3)*(-getindex(u,4) - getindex(u,5) + (1 + ^(getindex(u,3),2))*getindex(u,7)) + getindex(u,2)*(-getindex(u,4) - (1 + ^(getindex(u,3),2))*getindex(u,5) + ^(1 + ^(getindex(u,3),2),2)*getindex(u,7)))/(tau*(1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,(dtp*getindex(u,1)*getindex(u,3)*(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + getindex(u,2)*getindex(u,6) + getindex(u,3)*(getindex(u,4) + (1 + ^(getindex(u,2),2) + ^(getindex(u,3),2))*getindex(u,7)))/(tau*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,(-2*visc*(1 + ^(getindex(u,3),2))*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))/(3. *tau) + getindex(u,4)

,(4*visc*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)) + 3*tau*getindex(u,5))/(3. *^(tau,3))

,(-2*visc*getindex(u,2)*getindex(u,3)*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))/(3. *tau) + getindex(u,6)

,(zeta*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))/tau + getindex(u,7)

,(n + n*^(getindex(u,2),2) + getindex(u,2)*getindex(u,9) + getindex(u,3)*(n*getindex(u,3) + getindex(u,10)))/(tau*sqrt(1 + ^(getindex(u,2),2) + ^(getindex(u,3),2)))

,getindex(u,9)

,getindex(u,10)
)
return (At,Ax,Ay,source)
end


