# ==============================================================================
# Not usable: the second-moment rows 3-5 (pi_Q^r, pi_Q^perp, Pi_Q) of
# `one_d_viscous_matrix_fugacity_BG_only_second_moment` do not follow from the moment
# equation. Use `matrix1d_visc_HQ_BG_consistent_m2_5f!` (HQ_const_BG_consistent_m2.jl).
#
#   (1) The sigma_(nu) drive has the opposite sign: every nu-derivative coefficient is
#       minus the covariant one.
#   (2) The rows carry coordinate-transport couplings among (pi_l, pi_phi, Pi_Q) that are
#       absent from the Delta-projected D pi in the parallel-transported basis.
#   (3) Entries such as tauM*(2+tau^4)/3 and etaM*(-6+r^4+tau^4)/3 add tau^4 and r^4 to pure
#       numbers, which is dimensionally inconsistent. The moment block of A_t equals
#       tauM*u^tau*Identity only at tau = r = 1.
#
# Rows 1-2 (alpha, nu^r) are the grad-alpha-only first moment and are not affected.
#
# The code is kept to reproduce earlier results. The driver refuses to run unless forced with
# `allow_retracted = true` (closure keyword) or FLUIDUM_ALLOW_RETRACTED=1.
# `one_d_viscous_matrix_fugacity_BG_only_second_moment` itself stays callable: the consistent
# 5-field driver and the :legacy mode of the consistent second moment use it.
#
# The background (T, u) is prescribed, not evolved; `background_fields(tau, r)` returns
# (T, ur, dtT, drT, drur, dtur).
# ==============================================================================

function _refuse_retracted_second_moment(allow::Bool)
    allow && return nothing
    error("matrix1d_visc_HQ_BG_second_moment! is RETRACTED: its second-moment rows 3-5 " *
          "carry a flipped sigma_(nu) drive sign, spurious coordinate-transport couplings and " *
          "dimensionally inconsistent tau^4/r^4 entries (header of Matrix/HQ_const_BG_2nd_moment.jl). " *
          "Use matrix1d_visc_HQ_BG_consistent_m2_5f!. " *
          "To REPRODUCE an old artifact only, pass allow_retracted = true or set FLUIDUM_ALLOW_RETRACTED=1.")
end

function matrix1d_visc_HQ_BG_second_moment!(;dmn_eps=1e-6, background_fields, α_max=200.0, use_NR_tauM::Bool=false,
                                            allow_retracted::Bool = ALLOW_RETRACTED[])
    _refuse_retracted_second_moment(allow_retracted)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_second_moment!(A_i, Source, ϕ, tau, X, params;dmn_eps=dmn_eps,background_fields=background_fields, α_max=α_max, use_NR_tauM=use_NR_tauM, allow_retracted=true)
end


function matrix1d_visc_HQ_BG_second_moment!(A_i,Source,ϕ,tau,X,params;dmn_eps=1e-6,background_fields= nothing, α_max=200.0, use_NR_tauM::Bool=false,
                                            allow_retracted::Bool = ALLOW_RETRACTED[])
    _refuse_retracted_second_moment(allow_retracted)

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

    # τ_diffusion_hadron is bare: τ_n = D_s·I₃₁/(T·P₀), Eq. (30) of Capellino et al. 2205.07692,
    # degeneracy-free. τ_M and η_M are tied to `taun` below.
    # HQ_TAUN_SCALE is a diagnostic multiplier, default 1.0 (bare).
    taun=τ_diffusion_hadron(T,α_safe,params.eos,params.diffusion) * HQ_TAUN_SCALE[] #bare τ_n (Eq.30)
    #@show tau
    z = mq / T

    # Fokker-Planck to hydro matching:
    #   τ_M = (D_s/2)·[6z K1 + (z²+24) K2]/[z K1 + 4 K2] = (D_s z/2) K4/K3
    #   η_M = (D_s z/2) K3/K2 = (D_s/2)(4 + z K1/K2)   (the code carries an extra T by its σ/η convention,
    #         consistent across both branches: etaM_rel → T·D_s z/2 = T·taun/2 in the NR limit z→∞).
    # τ_n and τ_M are defined in the same convention, so their ratio is what is physical:
    #   τ_M/τ_n = ½·K4K2/K3² (≈ 0.550 at z = m_c/T_fo = 9.615).
    # Both are therefore tied to the bare `taun` above:
    #   τ_M = taun·½·K4K2/K3²,   η_M = T·taun/2  (the NR limit of the same tie).
    # The `_rel` forms below rebuild τ_M from the standalone `Ds` instead, which is consistent only
    # with a τ_n built from that same `Ds`; against the bare τ_n they give τ_M/τ_n ≈ 14.1 instead of
    # 0.550. They are a diagnostic branch, selected by env FIVO_IS2_TAUM_BARE=1.
    tauM_rel = Ds / (2) * (6 *z *besselk(1, z) + (z^2 + 24) * besselk(2, z))/(z*besselk(1,z) + 4 *besselk(2, z))
    etaM_rel = (Ds*T/2) * (4 + z *besselk(1, z)/besselk(2, z))
    tauM_deg = taun * besselk(4, z) * besselk(2, z) / (2 * besselk(3, z)^2)
    etaM_deg = T * taun / 2
    _tauM_bare = HQ_TAUM_BARE_DIAG[]          # env FIVO_IS2_TAUM_BARE, read once in __init__

    tauM = use_NR_tauM ? taun / 2     : (_tauM_bare ? tauM_rel : tauM_deg)
    etaM = use_NR_tauM ? T * taun / 2 : (_tauM_bare ? etaM_rel : etaM_deg)

    # FLUIDUM_HQ_MOMENTUM_DIMS=2: the 2-D-momentum-measure diagnostic, the same arithmetic as
    # `hq_transport` (EquationOfState/heavyquark.jl).
    if HQ_MOMENTUM_DIMS[] == 2
        λ1_2d = z * (z + 1) / (z^2 + 3z + 3)
        λ2_2d = z * (z^2 + 3z + 3) / (z^3 + 6z^2 + 15z + 15)
        taun *= (besselk(2, z) / besselk(3, z)) / λ1_2d
        tauM *= (besselk(3, z) / besselk(4, z)) / λ2_2d
        etaM  = T * taun / 2
    end

    # c_M: 0 by default (HQ_CM_BACKREACTION unset or "0"); D_s/T with HQ_CM_BACKREACTION=1.
    # Read once in __init__.
    cM   = HQ_CM_BACKREACTION[] ? Ds / T : zero(Ds)

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

    