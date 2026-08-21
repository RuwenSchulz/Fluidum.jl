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

    # τ_diffusion_hadron is BARE as of 2026-07-21 (its numerator now carries `Degeneracy`, which
    # cancels the one in `normalization`) — i.e. τ_n = D_s·I₃₁/(T·P₀), Eq. (30) of Capellino et al.
    # 2205.07692, degeneracy-free and equal to LangevInMedium.tau_n_main3. τ_M and η_M are tied to
    # `taun` below, so they follow automatically.
    # HQ_TAUN_SCALE is a diagnostic multiplier, DEFAULT 1.0 = bare (production). Set it to 1/6 only to
    # reproduce the historical ÷g_hq value for the convention study (diag_taun_hydro_vs_langevin.jl);
    # no production path sets it.
    taun=τ_diffusion_hadron(T,α_safe,params.eos,params.diffusion) * HQ_TAUN_SCALE[] #bare τ_n (Eq.30)
    #@show tau
    z = mq / T

    # Rigorous FP→hydro matching (Tex/HydroFPderivation/FP_Hydro_matching.tex, Eq. eq:tauM_etaM_expanded):
    #   τ_M = (D_s/2)·[6z K1 + (z²+24) K2]/[z K1 + 4 K2] = (D_s z/2) K4/K3
    #   η_M = (D_s z/2) K3/K2 = (D_s/2)(4 + z K1/K2)   (the code carries an extra T by its σ/η convention,
    #         consistent across both branches: etaM_rel → T·D_s z/2 = T·taun/2 in the NR limit z→∞).
    # The derivation defines τ_n AND τ_M in the SAME convention, so the RATIO is what is physical:
    #   τ_M/τ_n = ½·K4K2/K3² (≈ 0.550 at z = m_c/T_fo = 9.615).
    # PRODUCTION therefore ties both to `taun` — the bare Eq.(30) τ_n of line 29 — so the convention is
    # automatic and cannot drift when τ_n's normalisation is audited again:
    #   τ_M = taun·½·K4K2/K3²  (exact),   η_M = T·taun/2  (exact NR limit of the same tie).
    # The `_rel` forms below rebuild τ_M from the STANDALONE `Ds` of line 19 instead, which is only
    # consistent with a τ_n built from that same `Ds`; against the actual bare τ_n they give
    # τ_M/τ_n ≈ 14.1 rather than 0.550, i.e. 25.7× the tied value (4.28× after the historical ÷g_hq).
    # They are a diagnostic branch, NOT an alternative convention. Env FIVO_IS2_TAUM_BARE=1 selects them.
    # (An earlier version of this comment claimed τ_diffusion_hadron is degeneracy-weighted and that
    # τ_M/η_M must carry a matching ÷g_hq. That is false since 2026-07-21 — see line 22 — and the ÷g_hq
    # form is not what either branch computes. Do not reinstate it.)
    tauM_rel = Ds / (2) * (6 *z *besselk(1, z) + (z^2 + 24) * besselk(2, z))/(z*besselk(1,z) + 4 *besselk(2, z))
    etaM_rel = (Ds*T/2) * (4 + z *besselk(1, z)/besselk(2, z))
    tauM_deg = taun * besselk(4, z) * besselk(2, z) / (2 * besselk(3, z)^2)
    etaM_deg = T * taun / 2
    _tauM_bare = get(ENV, "FIVO_IS2_TAUM_BARE", "0") == "1"

    tauM = use_NR_tauM ? taun / 2     : (_tauM_bare ? tauM_rel : tauM_deg)
    etaM = use_NR_tauM ? T * taun / 2 : (_tauM_bare ? etaM_rel : etaM_deg)

    # ── FLUIDUM_HQ_MOMENTUM_DIMS=2 (2026-08-21): match the charm sector in the TWO-dimensional momentum
    # measure. The Langevin ensembles of LangevinPaper1 evolve transverse momenta only, so their
    # equilibrium is the 2-D Jüttner p dp e^{-E/T} and their ℓ-th moment relaxes at λ_ℓ(2D)·η_D, not at
    # the 3-D λ_ℓ = K_{ℓ+1}/K_{ℓ+2} that τ_n = D_s z K₃/K₂ and τ_M = D_s z K₄/(2K₃) encode. In that
    # measure every moment ∫_M^∞ E^n e^{-E/T} dE is an incomplete Gamma of integer order, so the ratios
    # are closed forms (checked against quadrature to 1e-9 at z = 2…20):
    #     λ₁(2D) = z(z+1)/(z²+3z+3),      λ₂(2D) = z(z²+3z+3)/(z³+6z²+15z+15),
    # and τ_n^{2D} = τ_n·(K₂/K₃)/λ₁(2D)  (5–12 % shorter over z = 3.5–10),  τ_M^{2D} = τ_M·(K₃/K₄)/λ₂(2D).
    # η_M stays tied to τ_n (η_M = T τ_n/2); κ = D_s n is dimension-independent. The 5/3 shear/bulk
    # tie and the tensor structure of the matrix are NOT switched (they are the 3-spatial-dimension
    # projection the solver is written in). Default "3" ⇒ production, bit-identical. DIAGNOSTIC knob:
    # it is a process-wide env read per call, so set it before the solve, never mid-run.
    if HQ_MOMENTUM_DIMS[] == 2
        λ1_2d = z * (z + 1) / (z^2 + 3z + 3)
        λ2_2d = z * (z^2 + 3z + 3) / (z^3 + 6z^2 + 15z + 15)
        taun *= (besselk(2, z) / besselk(3, z)) / λ1_2d
        tauM *= (besselk(3, z) / besselk(4, z)) / λ2_2d
        etaM  = T * taun / 2
    end

    # c_M = D_s/T is the rigorous 1st↔2nd-moment back-coupling (was hard-set to 0). Env-selectable so both
    # the no-backreaction (c_M=0) and the physical (c_M=D_s/T) cases can be compared:
    #   HQ_CM_BACKREACTION=0 → c_M=0 ;  otherwise (default) → c_M=D_s/T.
    cM   = get(ENV, "HQ_CM_BACKREACTION", "1") == "0" ? zero(Ds) : Ds / T

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

    