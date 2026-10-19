# ==============================================================================
# Thermodynamically consistent charm first moment on a prescribed 1+1D
# cylindrical background.
#
# The background (T, u) is prescribed, not evolved. `background_fields(tau, r)`
# returns (T, ur, dtT, drT, drur, dtur), as in HQ_const_BG.jl.
#
# Fields phi = (alpha, nu^r). The equation is the Delta-projected covariant
# divergence of the charm energy-momentum tensor of the rank-1 Grad (boosted
# Juttner) state,
#
#     (D_s/T) * Delta^r_nu Cov_mu T_Q^{mu nu} + nu^r = 0,
#     T_Q = epsQ u u + P0 Delta + h (u nu + nu u),   P0 = n T,  epsQ + P0 = n h,
#
# without the homogeneous-fluid step of Capellino et al. 2205.07692 App. D.
# That step (a) writes d^l P0 = T n0 d^l alpha by a constant-T identity and
# (b) drops all gradients of the thermodynamic coefficients. Derived
# symbolically, the consistent row is the grad-alpha-only row
# (kappa -> D_s dn/dalpha) plus five terms:
#
#       (D_s/T)(n + T dn/dT) * Delta^{r mu} d_mu T    pressure gradient,
#                                                     dP0 = n dmu + s dT
#       tau_n * n * a^r                               inertia, (eps + P0) a
#       tau_n * (nu . grad u)^r                       shear, vorticity, expansion
#       nu^r * tau_n * (theta_grad + D ln h)          expansion and the
#                                                     coefficient inside D
#       nu^r * tau_n * (u^tau/tau + u^r/r)            geometric dilution of the
#                                                     current (d(tau nu) = 0 in
#                                                     Bjorken flow)
#
# The geometric term comes from the Christoffel symbols of the divergence.
# All five are source terms (background gradients times fields); A_t and A_x
# are those of the grad-alpha-only row with kappa -> D_s dn/dalpha.
#
# Conventions: source on the LHS, A_t d_tau phi + A_x d_r phi + src = 0, as
# everywhere in this directory; the driver multiplies by A_t^{-1}. tau_n is
# `hq_taun` (Eq. 30 of 2205.07692). The consistent equation has the same
# relaxation time by construction: tau_n = D_s h/T defines h.
#
# Scope: first moment only. P0 = nT assumes unit charm charge per species;
# this is not checked. The 5-field variant below carries the second-moment rows
# of HQ_const_BG_2nd_moment.jl, which are not correct; at c_M = 0 they do not
# feed rows 1-2. The consistent second moment is HQ_const_BG_consistent_m2.jl.
# ==============================================================================

"""
    hqc_enthalpy_per_particle(T, α, eos, diffusion)

h = tau_n * T / D_s in GeV, the enthalpy per particle that tau_n encodes
(tau_n = D_s h / T). tau_n is `hq_taun` with both diagnostic knobs applied, so
the relation holds under every knob. Returns
(h, tau_n, D_s); `hq_transport(...)` returns the same three plus tau_M, eta_M.
"""
function hqc_enthalpy_per_particle(T, α, eos, diffusion)
    taun = hq_taun(T, α, eos, diffusion)
    Ds = DsT(diffusion, T) / T / fmGeV
    return taun * T / Ds, taun, Ds
end

"""
    hqc_dh_dT(T, α, eos, diffusion)

dh/dT by a central difference through `hqc_enthalpy_per_particle`, step
max(1e-4, 1e-3 T). EoS-agnostic: the closed form m K3'/K2' does not hold for
multi-species lists. Every consistent charm matrix uses this h'.
"""
function hqc_dh_dT(T, α, eos, diffusion)
    δ = max(1e-4, 1e-3 * T)
    return (hqc_enthalpy_per_particle(T + δ, α, eos, diffusion)[1] -
            hqc_enthalpy_per_particle(T - δ, α, eos, diffusion)[1]) / (2δ)
end

"""
    hqc_matrices(X, tau, r, ur, T, dtT, drT, drur, dtur,
                 n, dn_dT, dn_dmu, taun, Ds, h, hp; consistent = true)

The (alpha, nu^r) system. X = (alpha, nu^r). Returns (At, Ax, source) in the
source-on-the-LHS convention. `consistent = false` gives the grad-alpha-only
row of the homogeneous-fluid limit, with kappa = D_s * n (used in tests).
"""
function hqc_matrices(X, tau, r, ur, T, dtT, drT, drur, dtur,
                      n, dn_dT, dn_dmu, taun, Ds, h, hp; consistent::Bool = true)
    uτ  = sqrt(1 + ur^2)
    ν   = X[2]
    # the alpha-drive coefficient. ⚠ `dn_dmu` arrives regulated (+dmn_eps, from the
    # drivers), so the regulator also enters kappa = D_s dn/dalpha. Taking kappa
    # before the regulator changes nu^r at the 1e-3 level even where n > 0.1 n_max.
    κc  = consistent ? Ds * dn_dmu : Ds * n

    At = SMatrix{2,2}(
        uτ * dn_dmu,        # (1,1)  d_tau alpha  in row 1
        κc * ur * uτ,       # (2,1)  d_tau alpha  in row 2
        ur / uτ,            # (1,2)  d_tau nu     in row 1
        uτ * taun)          # (2,2)  d_tau nu     in row 2

    Ax = SMatrix{2,2}(
        ur * dn_dmu,        # (1,1)  d_r alpha  in row 1
        κc * uτ^2,          # (2,1)  d_r alpha  in row 2
        1.0,                # (1,2)  d_r nu     in row 1
        ur * taun)          # (2,2)  d_r nu     in row 2

    # row 1: continuity (as in the grad-alpha-only system; the 5-field normalization)
    src1 = n * (drur + dtur * ur / uτ + ur / r + uτ / tau) +
           ν * (1 / r + (ur + ur^3 + dtur * tau) / (uτ^3 * tau)) +
           (drT * ur + dtT * uτ) * dn_dT

    # row 2: the Delta-projection part common to both closures …
    src2 = ν * (1 - taun * (dtur * ur / uτ + drur * ur^2 / uτ^2))

    if consistent
        # … plus the five derived terms:
        #   pressure-gradient T channel + inertial   (field-independent sources)
        src2 += (Ds / T) * (drur * h * n * ur + dtur * h * n * uτ +
                            dtT * n * ur * uτ + dtT * dn_dT * T * ur * uτ +
                            drT * (n + dn_dT * T) * uτ^2)
        #   nu.grad u + expansion + D ln h           (proportional to nu)
        src2 += (Ds / T) * ν * (2 * drur * h + drT * hp * ur +
                                (2 * dtur * h * ur + dtT * hp * uτ^2) / uτ)
        #   geometric dilution of the current
        src2 += taun * ν * (uτ / tau + ur / r)
    end

    return (At, Ax, SVector{2}(src1, src2))
end

function matrix1d_visc_HQ_BG_consistent!(; dmn_eps = 1e-6, background_fields, α_max = 200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_consistent!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields, α_max = α_max)
end

function matrix1d_visc_HQ_BG_consistent!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0)

    T, ur, dtT, drT, drur, dtur = background_fields(tau, X[1])

    α_safe = clamp(ϕ[1], -α_max, α_max)

    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure                       # the charm-sector density (see charm_HRG.jl)
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu += dmn_eps

    h, taun, Ds = hqc_enthalpy_per_particle(T, α_safe, params.eos, params.diffusion)
    hp = hqc_dh_dT(T, α_safe, params.eos, params.diffusion)

    (At, Ax, source) = hqc_matrices(ϕ, tau, X[1], ur, T, dtT, drT, drur, dtur,
                                    n, dn_dT, dn_dmu, taun, Ds, h, hp; consistent = true)

    Ainv = inv(At)
    A_mul_B!(A_i[1], Ainv, Ax)
    jgemvavx!(Source, Ainv, source)
end

# ==============================================================================
# 5-field consistent system with the c_M back-coupling.
#
# phi = (alpha, nu^r, piQr, piQperp, PiQ). Rows 3-5 are the second-moment rows
# of `one_d_viscous_matrix_fugacity_BG_only_second_moment` (not correct, see
# HQ_const_BG_2nd_moment.jl; kept to reproduce earlier results), called with
# cM = 0 so that their own coupling entries do not enter. Rows 1-2 are the
# consistent first moment (hqc_matrices) plus the c_M force, derived
# symbolically with abstract indices:
#
#   cM * [ ur*ut*d_tau(zeta) + (1+ur^2)*d_r(zeta) + (drur*ur + dtur*ut)*zeta
#          + (1+ur^2)*(piQr - piQperp)/r + ur*ut*(2 piQr + piQperp)/tau ] ,
#   zeta = piQr + PiQ.
#
# The last two terms are the hoop-stress divergence and the tau redshift.
#
# cM = 0 by default (HQ_CM_BACKREACTION unset or "0"): rows 3-5 evolve and
# nothing feeds back, so (alpha, nu) equal the 2-field consistent solve.
# HQ_CM_BACKREACTION=1 gives cM = D_s/T, for which the system is elliptic.
# ==============================================================================

"""
    hqc_cm_row2(X, tau, r, ur, dtur, drur, cM) -> (at, ax, src)

The c_M force's row-2 pieces: `at` is the coefficient of BOTH d_tau piQr and
d_tau PiQ, `ax` that of both d_r piQr and d_r PiQ (the force acts on
zeta = piQr + PiQ), and `src` everything with no field derivative. piQperp has
NO derivative entries (the coupling carries none).
"""
function hqc_cm_row2(X, tau, r, ur, dtur, drur, cM)
    uτ = sqrt(1 + ur^2)
    piQr, piQperp, PiQ = X[3], X[4], X[5]
    ζ = piQr + PiQ
    at = cM * ur * uτ            # coefficient of d_tau piQr AND d_tau PiQ
    ax = cM * uτ^2               # coefficient of d_r  piQr AND d_r  PiQ
    src = cM * ((drur * ur + dtur * uτ) * ζ +
                uτ^2 * (piQr - piQperp) / r +
                ur * uτ * (2 * piQr + piQperp) / tau)
    return at, ax, src
end

"""
    hqc_first_moment_block(ϕ, tau, r, params, bg, dmn_eps, α_max, use_NR_tauM)

Everything the two 5-field consistent drivers (`matrix1d_visc_HQ_BG_consistent_5f!`
and `matrix1d_visc_HQ_BG_consistent_m2_5f!`) share, given the background
`bg = background_fields(tau, r)` = (T, ur, dtT, drT, drur, dtur): the charm thermodynamics (dn/dalpha regulated by `dmn_eps`), the transport
coefficients of `hq_transport`, h', rows 1-2 of the consistent first moment and
the c_M force's row-2 pieces. One body, so the two drivers' rows 1-2 stay identical.

c_M = D_s/T if `HQ_CM_BACKREACTION[]`, else 0 (the default), damped by
n/(n + `HQC_CM_NFLOOR[]`): the c_M force is not proportional to n, so in the
charm-free vacuum tail it drives nu >> n and collapses the implicit integrator's
dt on physics that carries no charm. The floor (default 1e-6) sits ~three orders
below n_fo ~ 8e-4, so at the freeze-out surface the suppression is < 0.2 %.
"""
@inline function hqc_first_moment_block(ϕ, tau, r, params, bg,
                                        dmn_eps, α_max, use_NR_tauM::Bool)
    # ⚠ the drivers evaluate `background_fields` THEMSELVES and pass the tuple: a
    # Function handed through a keyword-argument body without being called there is
    # not specialised on, and the dynamic call costs ~3.6 kB per cell.
    T, ur, dtT, drT, drur, dtur = bg
    α_safe = clamp(ϕ[1], -α_max, α_max)

    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu += dmn_eps

    tr = hq_transport(T, α_safe, params.eos, params.diffusion; use_NR_tauM = use_NR_tauM)
    hp = hqc_dh_dT(T, α_safe, params.eos, params.diffusion)

    At2, Ax2, src2 = hqc_matrices(ϕ, tau, r, ur, T, dtT, drT, drur, dtur,
                                  n, dn_dT, dn_dmu, tr.taun, tr.Ds, tr.h, hp; consistent = true)
    cM = HQ_CM_BACKREACTION[] ? tr.Ds / T : zero(tr.Ds)
    cM *= n / (n + HQC_CM_NFLOOR[])
    atc, axc, csrc = hqc_cm_row2(ϕ, tau, r, ur, dtur, drur, cM)
    return (T = T, ur = ur, dtT = dtT, drT = drT, drur = drur, dtur = dtur,
            n = n, dn_dT = dn_dT, dn_dmu = dn_dmu, tr = tr, hp = hp,
            At2 = At2, Ax2 = Ax2, src2 = src2, atc = atc, axc = axc, csrc = csrc)
end

function matrix1d_visc_HQ_BG_consistent_5f!(; dmn_eps = 1e-6, background_fields,
                                            α_max = 200.0, use_NR_tauM::Bool = false)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_consistent_5f!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields,
        α_max = α_max, use_NR_tauM = use_NR_tauM)
end

function matrix1d_visc_HQ_BG_consistent_5f!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0,
        use_NR_tauM::Bool = false)

    b = hqc_first_moment_block(ϕ, tau, X[1], params, background_fields(tau, X[1]),
                               dmn_eps, α_max, use_NR_tauM)

    # rows 3-5: the old second-moment rows (not correct, kept for reproducibility),
    # with their own cM set to zero; the c_M coupling enters through row 2 below
    Ats, Axs, srcs = one_d_viscous_matrix_fugacity_BG_only_second_moment(
        ϕ, tau, X[1], b.ur, b.T, b.dtT, b.drT, b.drur, b.dtur,
        b.n, b.dn_dmu, b.dn_dT, b.tr.taun, b.tr.Ds * b.n, b.tr.tauM, b.tr.etaM,
        params.diffusion.mass, zero(b.tr.Ds))

    # rows 1-2: consistent first moment + the c_M force
    At2, Ax2, src2 = b.At2, b.Ax2, b.src2
    At = MMatrix{5,5,Float64}(Ats); Ax = MMatrix{5,5,Float64}(Axs)
    src = MVector{5,Float64}(srcs)
    At[1, 1] = At2[1, 1]; At[1, 2] = At2[1, 2]; At[1, 3] = 0.0; At[1, 4] = 0.0; At[1, 5] = 0.0
    At[2, 1] = At2[2, 1]; At[2, 2] = At2[2, 2]; At[2, 3] = b.atc; At[2, 4] = 0.0; At[2, 5] = b.atc
    Ax[1, 1] = Ax2[1, 1]; Ax[1, 2] = Ax2[1, 2]; Ax[1, 3] = 0.0; Ax[1, 4] = 0.0; Ax[1, 5] = 0.0
    Ax[2, 1] = Ax2[2, 1]; Ax[2, 2] = Ax2[2, 2]; Ax[2, 3] = b.axc; Ax[2, 4] = 0.0; Ax[2, 5] = b.axc
    src[1] = src2[1]
    src[2] = src2[2] + b.csrc

    Ainv = inv(SMatrix{5,5}(At))
    A_mul_B!(A_i[1], Ainv, SMatrix{5,5}(Ax))
    jgemvavx!(Source, Ainv, SVector{5}(src))
end
