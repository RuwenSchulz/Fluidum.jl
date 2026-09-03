# the convention here are T, ux, uy, \[Pi]yy, \[Pi]zz, \[Pi]xy, \[Pi]B
#
# ─────────────────────────────────────────────────────────────────────────────────────────────
# 🔴 2026-09-03 — TEN ENTRIES OF THIS MATRIX ARE WRONG, AND THEY ARE ALL IN THE TWO SHEAR ROWS.
#
# Julia/tools/derive_2p1d_viscous.wls re-derives this system from scratch (plain Israel–Stewart,
# no second-order coefficients) and emits Matrix/viscous_generated.jl.  The SAME generator
# reproduces the shipped 1-D cylindrical `one_d_viscous_matrix` — the kernel every production
# 1-D run has used for a year — ENTRY BY ENTRY to 7.8e-16 (gate G4).  Against that reference the
# 7-field matrix below differs in exactly ten entries, over any state lattice:
#
#     A_t[4,3]  A_t[6,2]  A_t[6,3]    A_x[4,3]  A_x[6,2]  A_x[6,3]
#     A_y[4,2]  A_y[4,3]  A_y[6,2]  A_y[6,3]
#
# i.e. the velocity columns (2 = u^x, 3 = u^y) of the π^{yy} and π^{xy} rows.  Rows 1–3 (the
# conservation rows), row 5 (π^η_η), row 7 (Π) and the WHOLE source vector are correct, which is
# why Bjorken, the sound channel and the ideal cone all check out exactly.
#
# TWO SEPARATE DEFECTS LIVE IN THOSE TEN ENTRIES.
#  1. The σ^{xy} symmetrisation is lost.  At rest the row reads
#         τ_π D π^{xy} + π^{xy} = −2η ∂_x u^y        (A_y[6,2] ≡ 0 at every u, every π)
#     where Israel–Stewart asks for −2η σ^{xy} = −η(∂_x u^y + ∂_y u^x).  The SUM is right, the
#     SPLIT is not: for IRROTATIONAL transverse flow the two forms agree exactly, which is why
#     every smooth fireball started from rest — the FiVoHydro comparison included — has looked
#     right.  Under transverse VORTICITY it is wrong by −2η ω_{xy}: a rigid rotation, which must
#     produce no shear at all, drives π^{xy} at the full Navier–Stokes rate.  It also breaks
#     isotropy: at rest the transverse shear channel propagates at √(2 C_s) along x and does not
#     propagate at all along y, where both must be √(C_s).
#  2. Spurious τ_π·π velocity-gradient terms in rows 4 and 6.  At rest with η = 0 the derivation
#     gives zero for A_x[4,3], A_y[4,2], A_x[6,3], A_y[6,2]; this matrix gives ±τ_π π^{xy} and
#     ∓2τ_π π^{xy}.  It is NOT the co-rotating (Jaumann) coupling 2τ_π π_λ^{<μ}ω^{ν>λ}: adding
#     that term with EITHER sign makes the disagreement worse (12 entries instead of 10), and the
#     residual survives for irrotational gradients.  It is this defect that a row-by-row repair
#     cannot reach.
#
# Neither defect can be seen in 1-D: the transverse vorticity of a purely radial flow vanishes
# identically, so gate G4 passes either way.
#
# THE REPAIR IS OPT-IN AND DEFAULT OFF.  `FLUIDUM_2D_DERIVED=1`, or `Fluidum.VISC_2D_DERIVED[] =
# true`, routes `matrix2d_visc!` (and, through it, the 10-field `matrix2d_visc_HQ_BG!`, whose
# hydro block is this one) to `two_d_viscous_matrix_derived`.  With the flag off the path is
# bit-identical to every result produced before today.  Turning it on is RECOMMENDED for any new
# 2+1D work; it is left off so that nothing already on disk silently changes meaning.
# Gates: Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl §G, §R.
# ─────────────────────────────────────────────────────────────────────────────────────────────

"""
    VISC_2D_DERIVED

`true` when the 2+1D viscous kernels use `two_d_viscous_matrix_derived` — the matrix re-derived
from scratch by `Julia/tools/derive_2p1d_viscous.wls` — instead of the shipped
`two_d_viscous_matrix`, ten of whose entries are wrong (see the note at the top of this file).
Set from `FLUIDUM_2D_DERIVED` at load time; assignable at run time.
"""
const VISC_2D_DERIVED = Ref(get(ENV, "FLUIDUM_2D_DERIVED", "0") == "1")

"""
    two_d_viscous_matrix_active(u, tau, p, dtp, dtdtp, zeta, visc, tauS, tauB)

`two_d_viscous_matrix_derived` when `VISC_2D_DERIVED[]`, else the shipped
`two_d_viscous_matrix`.  Both have the same signature and the same field layout, so this is the
one place the choice is made.
"""
@inline two_d_viscous_matrix_active(args...) =
    VISC_2D_DERIVED[] ? two_d_viscous_matrix_derived(args...) : two_d_viscous_matrix(args...)

@inbounds @fastmath function matrix2d_visc!(A_i,Source,ϕ,t,X,params)
 

    therm=thermodynamic(ϕ[1],params)
        
    etaVisc=viscosity(ϕ[1],therm,params.shear)
    
    tauS=τ_shear(ϕ[1],therm,params.shear)
    
    tauB=τ_bulk(ϕ[1],therm,params.bulk)
    
    zeta=bulk_viscosity(ϕ[1],therm,params.bulk)
    
 
    (At,Ax, Ay, source)=two_d_viscous_matrix_active(ϕ,t,therm.pressure,therm.pressure_derivative[1],therm.pressure_hessian[1],zeta,etaVisc,tauS,tauB)
    

        
        Ainv= inv(At)

    
        A_mul_B!(A_i[1], Ainv,Ax)
        A_mul_B!(A_i[2], Ainv,Ay)
        
        jgemvavx!(Source, Ainv,source)
    
        
    
    end 



@inbounds @fastmath function two_d_viscous_matrix(u,tau,p,dtp,dtdtp,zeta,visc,tauS,tauB)
    
At=SMatrix{7,7}(

    dtp*(^(u[2],2) + ^(u[3],2)) + dtdtp*u[1]*(1 + ^(u[2],2) + ^(u[3],2))

    ,(dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(dtp + dtdtp*u[1])*u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,0

    ,0

    ,(2*(dtp*u[1]*u[2]*(1 + ^(u[3],2)) + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5])))/(1 + ^(u[3],2))

    ,(dtp*u[1]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2)) + (1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[7] + 3*u[2]*u[3]*u[6] + 2*^(u[2],3)*u[3]*u[6] + 3*u[2]*^(u[3],3)*u[6] - u[4] - 3*^(u[2],2)*u[4] - 2*^(u[2],4)*u[4] - ^(u[3],2)*u[4] - 3*^(u[2],2)*^(u[3],2)*u[4] - (1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[5])/((1 + ^(u[3],2))*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(u[6] + u[3]*(dtp*u[1]*u[2]*(1 + ^(u[2],2) + ^(u[3],2)) + ^(u[2],3)*u[7] + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4])))/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,(-2*(1 + ^(u[3],2))*(u[2]*visc + 3*tauS*u[3]*u[6]) + 6*tauS*u[2]*^(u[3],2)*u[4])/(3*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-2*u[2]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-3*u[2]*(-1 + tauS*(2 + 4*^(u[3],2)))*u[6] + 3*(1 + 2*tauS)*u[3]*u[4] + 6*tauS*(u[3] + ^(u[3],3))*u[5] + 2*^(u[2],2)*u[3]*(-2*visc + 3*tauS*(2*u[4] + u[5])))/(6*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(u[2]*zeta)/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(2*(u[2]*u[6] + u[3]*(dtp*u[1]*^(1 + ^(u[3],2),2) + ^(1 + ^(u[3],2),2)*u[7] + u[4] + u[2]*(-(u[3]*u[6]) + u[2]*(u[4] + u[5])))))/^(1 + ^(u[3],2),2)

    ,(dtp*u[1]*u[2]*u[3]*^(1 + ^(u[3],2),2)*(1 + ^(u[2],2) + ^(u[3],2)) - 2*^(u[2],4)*(-1 + ^(u[3],2))*u[6] + ^(1 + ^(u[3],2),2)*u[6] - 3*^(u[2],2)*(-1 + ^(u[3],4))*u[6] + 2*^(u[2],5)*u[3]*(u[4] + u[5]) + u[2]*u[3]*(1 + ^(u[3],2))*(^(1 + ^(u[3],2),2)*u[7] + 3*u[4] + u[5] + ^(u[3],2)*u[5]) + ^(u[2],3)*u[3]*(^(1 + ^(u[3],2),2)*u[7] + 5*u[4] + 3*u[5] + 3*^(u[3],2)*(u[4] + u[5])))/(^(1 + ^(u[3],2),2)*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(dtp*u[1]*(1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2)) + (1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2))*u[7] - u[2]*u[3]*u[6] + u[4] + ^(u[2],2)*u[4])/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,(4*^(u[3],3)*visc + 3*u[2]*u[6] + 6*tauS*u[2]*^(u[3],2)*u[6] + u[3]*((4 + 6*^(u[2],2))*visc + 3*u[4] - 6*tauS*(1 + ^(u[2],2))*u[4]))/(3*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-2*u[3]*visc)/(3*^(tau,2)*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-3*(-1 + 2*tauS)*u[3]*(1 + ^(u[3],2))*u[6] + 6*^(u[2],2)*u[3]*(1 + 2*tauS*^(u[3],2))*u[6] - 3*^(u[2],3)*(-4*(1 + ^(u[3],2))*visc + u[4] + 2*tauS*u[4] + u[5] + 2*tauS*^(u[3],2)*(2*u[4] + u[5])) + u[2]*(4*(3 + 5*^(u[3],2) + 2*^(u[3],4))*visc - 3*(u[4] + 2*tauS*u[4] + u[5]) - 3*^(u[3],2)*(u[5] + 2*tauS*(2*u[4] + u[5] + ^(u[3],2)*u[5]))))/(6*(1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(u[3]*zeta)/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(-^(u[2],2) + ^(u[3],2))/(1 + ^(u[3],2))

    ,-((u[2]*(1 + ^(u[2],2)))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2))))

    ,u[3]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,tauS*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,0

    ,-(^(u[2],2)/(1 + ^(u[3],2)))

    ,-((u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2)))/(1 + ^(u[3],2)))

    ,0

    ,0

    ,(tauS*sqrt(1 + ^(u[2],2) + ^(u[3],2)))/^(tau,2)

    ,0

    ,0

    ,(2*u[2]*u[3])/(1 + ^(u[3],2))

    ,(u[3] + 2*^(u[2],2)*u[3] + ^(u[3],3))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,u[2]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,tauS*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,^(u[2],2) + ^(u[3],2)

    ,u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,0

    ,0

    ,tauB*sqrt(1 + ^(u[2],2) + ^(u[3],2))

)
    #########################################################################################################################################################################
    
Ax=SMatrix{7,7}(

    (dtp + dtdtp*u[1])*u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,dtp + (dtp + dtdtp*u[1])*^(u[2],2)

    ,(dtp + dtdtp*u[1])*u[2]*u[3]

    ,0

    ,0

    ,0

    ,0

    ,(dtp*u[1]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2)) + (1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[7] + 3*u[2]*u[3]*u[6] + 2*^(u[2],3)*u[3]*u[6] + 3*u[2]*^(u[3],3)*u[6] - u[4] - 3*^(u[2],2)*u[4] - 2*^(u[2],4)*u[4] - ^(u[3],2)*u[4] - 3*^(u[2],2)*^(u[3],2)*u[4] - (1 + ^(u[2],2) + ^(u[3],2))*(1 + 2*^(u[2],2) + ^(u[3],2))*u[5])/((1 + ^(u[3],2))*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(2*(dtp*u[1]*u[2]*(1 + ^(u[3],2)) + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5])))/(1 + ^(u[3],2))

    ,u[3]*(dtp*u[1] + u[7])

    ,(-2*(1 + ^(u[3],2))*visc)/3. + (2*tauS*u[2]*u[3]*(-((1 + ^(u[3],2))*u[6]) + u[2]*u[3]*u[4]))/(1 + ^(u[2],2) + ^(u[3],2))

    ,(-2*visc)/(3*^(tau,2))

    ,(u[2]*(3*(1 - 2*tauS)*u[2]*u[6] - 12*tauS*u[2]*^(u[3],2)*u[6] + ^(u[3],3)*(-4*visc + 6*tauS*u[5]) + u[3]*(-4*(1 + ^(u[2],2))*visc + 3*(1 + 2*tauS)*u[4] + 6*tauS*(u[5] + ^(u[2],2)*(2*u[4] + u[5])))))/(6*(1 + ^(u[2],2) + ^(u[3],2)))

    ,zeta

    ,(dtp*u[1]*u[2]*u[3]*^(1 + ^(u[3],2),2)*(1 + ^(u[2],2) + ^(u[3],2)) - 2*^(u[2],4)*(-1 + ^(u[3],2))*u[6] + ^(1 + ^(u[3],2),2)*u[6] - 3*^(u[2],2)*(-1 + ^(u[3],4))*u[6] + 2*^(u[2],5)*u[3]*(u[4] + u[5]) + u[2]*u[3]*(1 + ^(u[3],2))*(^(1 + ^(u[3],2),2)*u[7] + 3*u[4] + u[5] + ^(u[3],2)*u[5]) + ^(u[2],3)*u[3]*(^(1 + ^(u[3],2),2)*u[7] + 5*u[4] + 3*u[5] + 3*^(u[3],2)*(u[4] + u[5])))/(^(1 + ^(u[3],2),2)*^(1 + ^(u[2],2) + ^(u[3],2),1.5))

    ,(2*(u[2]*(u[6] - ^(u[3],2)*u[6]) + u[3]*u[4] + ^(u[2],2)*u[3]*(u[4] + u[5])))/^(1 + ^(u[3],2),2)

    ,u[2]*(dtp*u[1] + u[7])

    ,u[6] + (2*u[2]*u[3]*((1 + ^(u[2],2) + ^(u[3],2))*visc + tauS*u[2]*u[3]*u[6]) - 2*tauS*(u[2] + ^(u[2],3))*u[3]*u[4])/(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,-(-4*(1 + ^(u[2],2))*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2))*visc + (-3 + 2*tauS)*u[2]*u[3]*(1 + ^(u[3],2))*u[6] - 2*^(u[2],3)*u[3]*(1 + 2*tauS*^(u[3],2))*u[6] + (1 + ^(u[3],2))*(2*u[4] + u[5] + ^(u[3],2)*u[5]) + ^(u[2],4)*(u[4] + tauS*(2 + 4*^(u[3],2))*u[4] + u[5] + 2*tauS*^(u[3],2)*u[5]) + ^(u[2],2)*((3 + 2*tauS + 2*(1 + 2*tauS)*^(u[3],2))*u[4] + 2*(1 + ^(u[3],2))*(1 + tauS*^(u[3],2))*u[5]))/(2*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)))

    ,0

    ,-((u[2]*(1 + ^(u[2],2)))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2))))

    ,-((1 + ^(u[2],2))/(1 + ^(u[3],2)))

    ,0

    ,tauS*u[2]

    ,0

    ,0

    ,0

    ,-((u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2)))/(1 + ^(u[3],2)))

    ,-1 - ^(u[2],2)/(1 + ^(u[3],2))

    ,0

    ,0

    ,(tauS*u[2])/^(tau,2)

    ,0

    ,0

    ,(u[3] + 2*^(u[2],2)*u[3] + ^(u[3],3))/((1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(2*u[2]*u[3])/(1 + ^(u[3],2))

    ,1

    ,0

    ,0

    ,tauS*u[2]

    ,0

    ,u[2]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,1 + ^(u[2],2)

    ,u[2]*u[3]

    ,0

    ,0

    ,0

    ,tauB*u[2]

)
    #########################################################################################################################################################################
Ay=SMatrix{7,7}(

    (dtp + dtdtp*u[1])*u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,(dtp + dtdtp*u[1])*u[2]*u[3]

    ,dtp + (dtp + dtdtp*u[1])*^(u[3],2)

    ,0

    ,0

    ,0

    ,0

    ,(u[6] + u[3]*(dtp*u[1]*u[2]*(1 + ^(u[2],2) + ^(u[3],2)) + ^(u[2],3)*u[7] + u[3]*u[6] + u[2]*(u[7] + ^(u[3],2)*u[7] - u[4])))/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,u[3]*(dtp*u[1] + u[7])

    ,0

    ,-(((1 + 2*tauS*^(u[3],2))*(u[6] + ^(u[3],2)*u[6] - u[2]*u[3]*u[4]))/(1 + ^(u[2],2) + ^(u[3],2)))

    ,0

    ,(-(u[2]*u[3]*(1 + tauS*(2 + 4*^(u[3],2)))*u[6]) + 2*u[4] + u[5] + ^(u[2],2)*(1 + 2*tauS*^(u[3],2))*(2*u[4] + u[5]) + ^(u[3],2)*(u[4] + 2*tauS*u[4] + u[5] + 2*tauS*(1 + ^(u[3],2))*u[5]))/(2*(1 + ^(u[2],2) + ^(u[3],2)))

    ,0

    ,(dtp*u[1]*(1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2)) + (1 + ^(u[2],2) + ^(u[3],2))*(1 + ^(u[2],2) + 2*^(u[3],2))*u[7] - u[2]*u[3]*u[6] + u[4] + ^(u[2],2)*u[4])/^(1 + ^(u[2],2) + ^(u[3],2),1.5)

    ,u[2]*(dtp*u[1] + u[7])

    ,2*u[3]*(dtp*u[1] + u[7])

    ,(4*(1 + ^(u[3],2))*visc)/3. + (u[3]*(u[2]*(1 + 2*tauS*^(u[3],2))*u[6] + (1 - 2*tauS*(1 + ^(u[2],2)))*u[3]*u[4]))/(1 + ^(u[2],2) + ^(u[3],2))

    ,(-2*visc)/(3*^(tau,2))

    ,(u[3]*(-3*(-1 + 2*tauS)*u[3]*(1 + ^(u[3],2))*u[6] + 6*^(u[2],2)*u[3]*(1 + 2*tauS*^(u[3],2))*u[6] + ^(u[2],3)*(8*(1 + ^(u[3],2))*visc - 6*tauS*^(u[3],2)*(2*u[4] + u[5]) - 3*(u[4] + 2*tauS*u[4] + u[5])) + u[2]*(8*^(1 + ^(u[3],2),2)*visc - 3*(u[4] + tauS*(2 + 4*^(u[3],2))*u[4] + (1 + ^(u[3],2))*(1 + 2*tauS*^(u[3],2))*u[5]))))/(6*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)))

    ,zeta

    ,u[3]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,0

    ,1

    ,tauS*u[3]

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,0

    ,(tauS*u[3])/^(tau,2)

    ,0

    ,0

    ,u[2]/sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,1

    ,0

    ,0

    ,0

    ,tauS*u[3]

    ,0

    ,u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))

    ,u[2]*u[3]

    ,1 + ^(u[3],2)

    ,0

    ,0

    ,0

    ,tauB*u[3]

)
    #########################################################################################################################################################


    #########################################################################################################################################################

source=SVector{7}(

    (dtp*u[1]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)) + u[7] + 2*u[2]*u[3]*u[6] + ^(u[2],2)*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5]) + u[5] + ^(u[3],2)*((2 + ^(u[3],2))*u[7] + u[4] + u[5]))/(tau*(1 + ^(u[3],2)))

    ,(dtp*u[1]*u[2]*(1 + ^(u[3],2))*(1 + ^(u[2],2) + ^(u[3],2)) + 2*^(u[2],2)*u[3]*u[6] + u[3]*(1 + ^(u[3],2))*u[6] + ^(u[2],3)*(u[7] + ^(u[3],2)*u[7] - u[4] - u[5]) + u[2]*(^(1 + ^(u[3],2),2)*u[7] - u[4] - (1 + ^(u[3],2))*u[5]))/(tau*(1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(dtp*u[1]*u[3]*(1 + ^(u[2],2) + ^(u[3],2)) + ^(u[3],3)*u[7] + u[2]*u[6] + u[3]*(u[7] + ^(u[2],2)*u[7] + u[4]))/(tau*sqrt(1 + ^(u[2],2) + ^(u[3],2)))

    ,(-2*(1 + ^(u[3],2))*sqrt(1 + ^(u[2],2) + ^(u[3],2))*visc)/(3*tau) + u[4]

    ,(4*sqrt(1 + ^(u[2],2) + ^(u[3],2))*visc + 3*tau*u[5])/(3*^(tau,3))

    ,(-2*u[2]*u[3]*sqrt(1 + ^(u[2],2) + ^(u[3],2))*visc)/(3*tau) + u[6]

    ,(sqrt(1 + ^(u[2],2) + ^(u[3],2))*zeta)/tau + u[7]

)

return (At,Ax,Ay,source)
end


