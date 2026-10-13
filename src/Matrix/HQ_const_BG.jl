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
#
# FIELDS: phi = (alpha, nu^r) -- two fields, the charm first moment only (the medium is the closure
# `background_fields`).  Until 2026-10-13 this line listed the 7-field layout of
# 1d_viscous_HQ_cilindrical_fugacity.jl, copied with the code.
# EQUATIONS: Nabla.(n u + nu) = 0;  tau_n Delta D nu + nu + kappa Delta grad(alpha) = 0 (the shipped,
# grad-alpha-only closure).  GATED against a covariant derivation, read out of this file: gate M-BG of
# Julia/tools/derive_1d_matrices.wls (raw entries, symbolically) and J-BG of
# Projects/FluidumValidation/gate_matrices_1d.jl (in programme.jl check); also H-G4 of
# derive_hq_2p1d.wls and G3 of Tex/MaxEntHydro/diag_hq_consistent_gates.jl.




function matrix1d_visc_HQ_BG!(;dmn_eps=1e-6, background_fields, α_max=200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields, α_max=α_max)
end


#THIS IS THE MATRIX THAT DID NOT CREATE PROBLEMS WITH THE BUMP (GUBSER)
# dmn_eps: 1e-6 in BOTH methods since 2026-10-13. The direct method defaulted to 1e-10 while the
# closure above (the only route production takes: every mode-registry and driver call builds the
# closure, and the closure always passes its own value through) defaulted to 1e-6 -- so 1e-10 was
# reachable only by calling this method by hand without the keyword. No production output changes.
function matrix1d_visc_HQ_BG!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-6,background_fields= nothing, α_max=200.0)

    T,ur,dtT,drT,drur,dtur = background_fields(tau,X[1])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    
    thermo = thermodynamic(T,α_safe,params.eos.hadron_list)
    n=thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu+= dmn_eps

    κ = diffusion_hadron(T,α_safe,params.eos,params.diffusion) #diffusion coefficient for hadrons
    # τ_n with the diagnostic knobs (HQ_TAUN_SCALE, FLUIDUM_HQ_MOMENTUM_DIMS) applied, as in every
    # charm matrix (hq_taun, 2026-10-13; this file ignored both before). Production: bare τ_n.
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


