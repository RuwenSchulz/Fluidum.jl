# ⚠ "const_BG" IS A MISNOMER, AND THE NAME IS KEPT ONLY BECAUSE IT IS ON DISK.
# The background here is NOT constant. `background_fields(tau, X[1])` returns
#     T, ur, dtT, drT, drur, dtur
# at that (tau, r) — a fully time- and space-dependent medium, splined over an r-grid AND a
# tau-grid (see Projects/LangevinPaper1/generate_background_splines.jl). What "BG" actually means
# is PRESCRIBED rather than CO-EVOLVED: the medium is supplied to the charm sector instead of being
# advanced alongside it.
#
# 🔑 AND THAT IS NOT AN APPROXIMATION. At c_M = 0 the charm is a PROBE — it does not back-react on
# the medium. Measured 2026-09-15 (Projects/FiVoFluidumComparison/COMPARISON_1P1D.md §3.2), two ways:
#   * structure: in the live 7-field matrix1d_visc_HQ!, max|B[1:5, 6:7]| = 0.000e+00 exactly, and no
#     medium source row references alpha or nu^r;
#   * behaviour: quadrupling D_s (0.1 -> 0.4) moves nu^r by 62 % and the medium by 6.6e-16 ... 8.6e-15.
# So prescribing the medium and co-evolving it are the SAME problem, up to discretisation. Solving
# it in two stages is a choice about scheduling, not about physics.
#
# ⛔ The symbols `:HQ_const_BG*` are written into stored artifacts as `evo_mode`; renaming them
# would orphan every .jld2 that carries one. The name stays; this note is the fix.
# the convention here are T, ur,  \[Pi]phiphi, \[Pi]etaeta, \[Pi]B, mu, nur




function matrix1d_visc_HQ_BG!(;dmn_eps=1e-6, background_fields, α_max=200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields, α_max=α_max)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
function matrix1d_visc_HQ_BG!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-10,background_fields= nothing, α_max=200.0)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    
    thermo = thermodynamic(T,α_safe,params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu+= dmn_eps

    κ = diffusion_hadron(T,α_safe,params.eos,params.diffusion) #diffusion coefficient for hadrons
    tauDiffusion=τ_diffusion_hadron(T,α_safe,params.eos,params.diffusion) #tau diffusion for hadrons

    #actually our equations don t depend on p: we can just put as entry dpt instead, in any case it will not be used (but in the future maybe it will be )
    #(At,Ax, source)=one_d_viscous_HQ_matrix(ϕ,t,X[1],dpt,dpt,dptt,zeta,etaVisc,tauS,tauB,n,dtn,dmn,tauDiff,Ds)
    (At,Ax, source)=one_d_viscous_matrix_fugacity_BG_only(ϕ,tau,X[1],ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,κ)

        
    Ainv= inv(At)


    A_mul_B!(A_i[1], Ainv,Ax)
    
    
    jgemvavx!(Source, Ainv,source)

    end 



function one_d_viscous_matrix_fugacity_BG_only(X,tau,r,ur,dtT,drT,drur,dtur,n,dn_dmu,dn_dT,tauDiffusion,κ)
    #these are the matrices from 1d viscous hydro
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


