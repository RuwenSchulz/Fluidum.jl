# ==============================================================================
# HQ_const_BG_consistent_m2.jl — THE THERMODYNAMICALLY CONSISTENT SECOND MOMENT.
#
# NEW FILE (2026-09-02). NOT included by src/Fluidum.jl — nothing shipped
# changes. Drivers load it into the module at runtime AFTER
# HQ_const_BG_consistent.jl (whose hqc_matrices / hqc_cm_row2 /
# hqc_enthalpy_per_particle this file calls):
#
#     Base.include(Fluidum, ".../Matrix/HQ_const_BG_consistent.jl")
#     Base.include(Fluidum, ".../Matrix/HQ_const_BG_consistent_m2.jl")
#
# and register :HQ_const_BG_consistent_m2_5f themselves.
#
# WHAT IT IS. Rows 1-2 (alpha, nu^r) are BYTE-IDENTICAL in structure to
# :HQ_const_BG_consistent_5f (same hqc_matrices call, same cM plumbing). Rows
# 3-5 replace the shipped nabla-nu-only second-moment relaxation with the
# thermodynamically consistent equations of
# Tex/LangevinPaper1/M2_CONSISTENT_DERIVATION.md (Sec. 4, delta-M-extension
# convention — the variant the kinetic validation favours 6-8x over strict-Pi),
# reduced under the paper's symmetry: boost-invariant, azimuthally symmetric
# Milne, u = (u^tau, u^r, 0, 0), mostly-plus g = diag(-1, 1, r^2, tau^2).
#
# FIELD CONVENTION (matches the verified c_M coupling of
# derive_hq_cm_coupling.wls and the Langevin extraction): the transverse basis
# is the PARALLEL-TRANSPORTED orthonormal triad
#     l    = (u^r, u^tau, 0, 0)          (unit, ⊥u, in the tau-r plane)
#     phî  = (0, 0, 1/r, 0),   etâ = (0, 0, 0, 1/tau)
# for which  D l = a_l u  (so Delta·Dl = 0)  and  D phî = D etâ = 0, hence
# Delta-projected D pi is DIAGONAL with entries D p_i.  Fields:
#     X3 = p_l   (pi_Q in the l̂ l̂ channel; = pi_Q^r_r at the axis)
#     X4 = p_phi (pi_Q^phi_phi),     p_eta = -(X3 + X4)  (traceless)
#     X5 = Pi_Q  (trace channel),    dM^i_i = p_i + Pi_Q
#
# THE REDUCED EQUATIONS (all-on-LHS = 0, the file convention
# At d_tau phi + Ax d_r phi + src = 0; D = u^tau d_tau + u^r d_r,
# D_l = u^r d_tau + u^tau d_r on scalars).  Kinematics:
#     theta_l = D_l u^r / u^tau,  theta_phi = u^r/r,  theta_eta = u^tau/tau
#     theta   = theta_l + theta_phi + theta_eta,   s_i = theta_i - theta/3
#     a_l     = (D u^r)/u^tau,    nu_l = nu^r/u^tau,   (vorticity ⊥-block = 0)
#     (grad nu)_ll   = D_l nu^r / u^tau - u^r nu^r (D_l u^r)/u^tau^3
#     (grad nu)_phph = nu^r/r,      (grad nu)_etet = u^r nu^r/(u^tau tau)
#     theta_nu = d_tau nu^tau + d_r nu^r + nu^tau/tau + nu^r/r,
#     a.nu = a_l nu_l;   identity (checked): sum_i (grad nu)_ii = theta_nu - a.nu
#     signu_i = (grad nu)_ii - (theta_nu - a.nu)/3        [sigma_(nu) diagonal]
#
# Traceless channels i in {l, phi}  (E_eta = -(E_l + E_phi) identically):
#     tau_M D p_i + p_i + 2 eta_Q signu_i                              (legacy)
#       + tau_M [ ((5/3)theta + D ln C) p_i + 2(p_i s_i - (pi:sigma)/3)
#                 + 2 Pi_Q s_i ]                                       (iii + dM-ext)
#       + 2 etabar s_i                                                 (ii)
#       + [ 2 lambda_a a_l nu_l + (D_s/T) nu_l D_l(Th) ] (delta_il - 1/3)  (iv)
#     = 0,          pi:sigma = p_l s_l + p_phi s_phi + p_eta s_eta
#
# Trace channel:
#     tau_M D Pi_Q + Pi_Q + zeta_Q theta_nu                            (legacy)
#       + tau_M ((5/3)theta + D ln C) Pi_Q + (2/3) tau_M pi:sigma      (iii + dM-ext)
#       + (tau_n P0/2) [ D alpha + (A/B) D ln T + (5/3) theta ]        (ii)
#       + (D_s A/3T) a.nu + (5 D_s/6T) nu_l D_l(Th)                    (iv)
#     = 0
#
# COEFFICIENTS (paper values; z = m_q/T, K_n = besselk):
#     tau_M  = tau_n K4 K2 / (2 K3^2)        eta_Q  = T tau_n / 2   (= etaM)
#     zeta_Q = (5/3) eta_Q                   etabar = tau_n P0 / 2,  P0 = n T
#     lambda_a = (D_s/2T)(m^2 + 6 T h)       A = m^2 + 5 T h,  A/B = 5 + z K2/K3
#     D ln C = [z ((K3+K5)/(2K4) - (K2+K4)/(2K3))] D ln T,   K5 = K3 + 8 K4/z
#     grad(Th) = (h + T h') grad T           (h EOS-agnostic per hqc precedent)
# For a multi-species hadron list, h and n are the list's aggregates (exactly as
# rows 1-2 and the shipped tau_M plumbing) while the m^2 terms use
# params.diffusion.mass — the same convention tau_M itself has always used.
#
# THE SOURCE SWITCH — env LP1_M2_SOURCES, read per call (like
# HQ_CM_BACKREACTION) so no const-capture can pin it at include time:
#
#   "consistent" (default)  the equations above, all completion terms on.
#   "legacy"                rows 3-5 delegated VERBATIM to the shipped
#                           one_d_viscous_matrix_fugacity_BG_only_second_moment
#                           (cM = 0), i.e. bit-identical to
#                           :HQ_const_BG_consistent_5f.  Gate G1's control: an
#                           LP1_M2_SOURCES=legacy solve must reproduce the
#                           charm_hydro_consistent payloads at solver roundoff.
#   "nabla_nu"              the consistent rows with the NEW terms zeroed —
#                           the paper's nabla-nu-only relaxation
#                           tau_M D pi + pi = -2 eta_Q sigma_(nu),
#                           tau_M D Pi + Pi = -zeta_Q theta_nu
#                           in the clean covariant reduction (= the derivation
#                           note's "shipped ODE" curves at the axis).
#
# 🔴 WHY "legacy" IS A DELEGATION AND NOT THE nabla_nu LIMIT (2026-09-02,
# measured by rotating the shipped rows into this file's row basis —
# scratch diag_shipped_row_content.jl):  the shipped Mathematica rows are NOT
# the nabla-nu-only limit of the covariant reduction.  Two structural
# differences, both exact:
#   (1) the sigma_(nu) drive enters the shipped rows with the OPPOSITE SIGN
#       (every nu-derivative coefficient is exactly minus the covariant one —
#       the mostly-minus notebook heritage, same family as the W0/C0 traps);
#   (2) the shipped rows carry coordinate-transport couplings among
#       (pi_l, pi_phi, Pi_Q) (e.g. d src3/d X3 = 0.842 tau_M-units instead of
#       1, cross terms ±0.16, 0.23) that the parallel-transported basis shows
#       to be absent from Delta-proj D pi.
# So "zeroing the new terms" cannot reproduce the production payloads; only the
# verbatim shipped rows can, and the consistent-vs-legacy difference contains
# the drive-sign repair and the transport cleanup ON TOP of classes ii-iv.
# nabla_nu exists precisely to attribute the split.
#
# c_M SEMANTICS UNCHANGED: rows 3-5 are PASSIVE (production runs
# HQ_CM_BACKREACTION=0 ⇒ cM = 0, and rows 3-5 carry no cM entries of their
# own); the row-2 coupling is hqc_cm_row2 exactly as in the first-moment file.
#
# Derivation + validation: Tex/LangevinPaper1/M2_CONSISTENT_DERIVATION.md
# (symbolic .wls identity exact; manufactured-field FD residual 1.5e-11;
# kinetic validation mean|dev| tau in [5,13] = 0.0087/0.0098 const/linear).
# The general-r reduction above is verified pointwise against the covariant
# E+N+P assembly in the scratch gate accompanying the 2026-09-02 session
# (see M2_CONSISTENT_DERIVATION.md App. A) and against the shipped rows via
# the At^-1 Ax / At^-1 src equivalence test in legacy mode.
# ==============================================================================

"""
    hqc_m2_sources() -> Symbol

The LP1_M2_SOURCES switch, read per call: :consistent (default), :legacy
(shipped rows verbatim — the G1 control) or :nabla_nu (consistent minus the
new terms). Anything else is an error — a misspelled control is a wrong
physics run.
"""
function hqc_m2_sources()
    s = get(ENV, "LP1_M2_SOURCES", "consistent")
    s == "consistent" && return :consistent
    s == "legacy" && return :legacy
    s == "nabla_nu" && return :nabla_nu
    error("LP1_M2_SOURCES must be \"consistent\", \"legacy\" or \"nabla_nu\", got \"$s\"")
end

"""
    hqc_m2_rows(X, tau, r, ur, T, dtT, drT, drur, dtur,
                n, taun, Ds, h, hp, tauM, etaM, mq; sources = :consistent)

Rows 3-5 of the consistent-second-moment system: returns
(at3, ax3, src3) as 3x5 SMatrix rows / length-3 source, in the
source-on-the-LHS convention. `sources = :nabla_nu` zeroes exactly the new
terms; `:legacy` is handled one level up (verbatim shipped rows — see the
file header).
"""
function hqc_m2_rows(X, tau, r, ur, T, dtT, drT, drur, dtur,
                     n, taun, Ds, h, hp, tauM, etaM, mq;
                     sources::Symbol = :consistent)
    uτ  = sqrt(1 + ur^2)
    ν   = X[2]
    p_l = X[3]; p_φ = X[4]; Π = X[5]
    p_η = -(p_l + p_φ)

    # ── kinematics of the background flow ─────────────────────────────────────
    Dlur = ur * dtur + uτ * drur            # D_l u^r
    θ_l  = Dlur / uτ
    θ_φ  = ur / r
    θ_η  = uτ / tau
    θ    = θ_l + θ_φ + θ_η
    s_l  = θ_l - θ / 3
    s_φ  = θ_φ - θ / 3
    s_η  = θ_η - θ / 3
    a_l  = (uτ * dtur + ur * drur) / uτ     # a·l = (D u^r)/u^tau
    ν_l  = ν / uτ

    # ── sigma_(nu): derivative coefficients (columns 2) and nu-proportional parts
    #     d_tau nu^r enters (grad nu)_ll and theta_nu with coefficient u^r/u^tau,
    #     d_r nu^r with coefficient 1  ⇒  signu_l carries 2/3 of each, signu_{phi,eta} -1/3.
    cν_t = ur / uτ                          # d_tau nu coefficient of (grad nu)_ll and theta_nu
    Sll  = -ur * ν * Dlur / uτ^3            # nu-proportional part of (grad nu)_ll
    Sφφ  = ν / r
    Sηη  = ur * ν / (uτ * tau)
    aν   = a_l * ν_l
    Strν_noa = ν * dtur / uτ^3 + ν / r + ur * ν / (uτ * tau)   # nu-prop. part of theta_nu
    Strν = Strν_noa - aν                                        # ... of (theta_nu - a·nu)

    ηQ = etaM                               # = T taun/2 — the paper's eta_Q
    ζQ = (5.0 / 3.0) * ηQ

    # ── legacy rows: tau_M D p_i + p_i = -2 eta_Q signu_i ; trace with zeta_Q theta_nu
    at3 = zeros(MMatrix{3,5,Float64})
    ax3 = zeros(MMatrix{3,5,Float64})
    src = zeros(MVector{3,Float64})

    # relaxation advection
    at3[1, 3] = tauM * uτ;  ax3[1, 3] = tauM * ur
    at3[2, 4] = tauM * uτ;  ax3[2, 4] = tauM * ur
    at3[3, 5] = tauM * uτ;  ax3[3, 5] = tauM * ur

    # sigma_(nu) drives (derivative entries in column 2)
    at3[1, 2] = 2 * ηQ * (2.0 / 3.0) * cν_t;  ax3[1, 2] = 2 * ηQ * (2.0 / 3.0)
    at3[2, 2] = -2 * ηQ * (1.0 / 3.0) * cν_t; ax3[2, 2] = -2 * ηQ * (1.0 / 3.0)
    at3[3, 2] = ζQ * cν_t;                    ax3[3, 2] = ζQ

    src[1] = p_l + 2 * ηQ * (Sll - Strν / 3)
    src[2] = p_φ + 2 * ηQ * (Sφφ - Strν / 3)
    src[3] = Π + ζQ * Strν_noa

    sources == :legacy && error("hqc_m2_rows: :legacy is handled by the caller (verbatim shipped rows)")
    if sources == :consistent
        z  = mq / T
        K2 = besselk(2, z); K3 = besselk(3, z); K4 = besselk(4, z)
        K5 = K3 + 8 * K4 / z
        dlnC_dlnT = z * ((K3 + K5) / (2 * K4) - (K2 + K4) / (2 * K3))
        DlnT = (uτ * dtT + ur * drT) / T
        DlnC = dlnC_dlnT * DlnT
        DlTh = (h + T * hp) * (ur * dtT + uτ * drT)     # D_l(Th)

        ηbar = taun * n * T / 2                          # tau_n P0 / 2
        λa   = Ds * (mq^2 + 6 * T * h) / (2 * T)         # (Ds/2T)(A+B)
        Aco  = mq^2 + 5 * T * h                          # A = I41/P0
        ABr  = 5 + z * K2 / K3                           # A/B = dln I31/dln T
        rank1 = 2 * λa * a_l * ν_l + (Ds / T) * ν_l * DlTh   # class-(iv) traceless bracket

        πσ = p_l * s_l + p_φ * s_φ + p_η * s_η
        geo = tauM * ((5.0 / 3.0) * θ + DlnC)

        # traceless channels: (iii) + dM-extension + (ii) + (iv)
        src[1] += geo * p_l + 2 * tauM * (p_l * s_l - πσ / 3) + 2 * tauM * Π * s_l +
                  2 * ηbar * s_l + rank1 * (2.0 / 3.0)
        src[2] += geo * p_φ + 2 * tauM * (p_φ * s_φ - πσ / 3) + 2 * tauM * Π * s_φ +
                  2 * ηbar * s_φ + rank1 * (-1.0 / 3.0)
        # trace channel: (iii) + dM-extension + (ii: background bulk incl. D alpha) + (iv)
        src[3] += geo * Π + (2.0 / 3.0) * tauM * πσ +
                  ηbar * (ABr * DlnT + (5.0 / 3.0) * θ) +
                  (Ds * Aco / (3 * T)) * aν + (5.0 / 6.0) * (Ds / T) * ν_l * DlTh
        # the D alpha of D ln I31 — matrix entries in column 1 of the trace row
        at3[3, 1] = ηbar * uτ
        ax3[3, 1] = ηbar * ur
    end

    return SMatrix{3,5}(at3), SMatrix{3,5}(ax3), SVector{3}(src)
end

function matrix1d_visc_HQ_BG_consistent_m2_5f!(; dmn_eps = 1e-6, background_fields,
                                               α_max = 200.0, use_NR_tauM::Bool = false)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_consistent_m2_5f!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields,
        α_max = α_max, use_NR_tauM = use_NR_tauM)
end

function matrix1d_visc_HQ_BG_consistent_m2_5f!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0,
        use_NR_tauM::Bool = false)

    T, ur, dtT, drT, drur, dtur = background_fields(tau, X[1])
    α_safe = clamp(ϕ[1], -α_max, α_max)

    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_dmu = thermo.pressure_derivative
    dn_dmu += dmn_eps

    h, taun, Ds = hqc_enthalpy_per_particle(T, α_safe, params.eos, params.diffusion)
    δ  = max(1e-4, 1e-3 * T)
    hp = (hqc_enthalpy_per_particle(T + δ, α_safe, params.eos, params.diffusion)[1] -
          hqc_enthalpy_per_particle(T - δ, α_safe, params.eos, params.diffusion)[1]) / (2δ)

    z = params.diffusion.mass / T
    tauM = use_NR_tauM ? taun / 2 : taun * besselk(4, z) * besselk(2, z) / (2 * besselk(3, z)^2)
    etaM = T * taun / 2

    # rows 3-5: the consistent second moment; under LP1_M2_SOURCES=legacy the
    # VERBATIM shipped rows (gate G1's control — see the header on why the
    # shipped rows are not the nabla_nu limit); under nabla_nu the clean
    # covariant nabla-nu-only forms.
    src_mode = hqc_m2_sources()
    local at3, ax3, src3
    if src_mode == :legacy
        kappa = Ds * n
        Ats, Axs, srcs = one_d_viscous_matrix_fugacity_BG_only_second_moment(
            ϕ, tau, X[1], ur, T, dtT, drT, drur, dtur,
            n, dn_dmu, dn_dT, taun, kappa, tauM, etaM, params.diffusion.mass, zero(Ds))
        at3 = SMatrix{3,5,Float64}(ntuple(k -> Ats[3 + (k - 1) % 3, 1 + (k - 1) ÷ 3], 15))
        ax3 = SMatrix{3,5,Float64}(ntuple(k -> Axs[3 + (k - 1) % 3, 1 + (k - 1) ÷ 3], 15))
        src3 = SVector{3,Float64}(srcs[3], srcs[4], srcs[5])
    else
        at3, ax3, src3 = hqc_m2_rows(ϕ, tau, X[1], ur, T, dtT, drT, drur, dtur,
                                     n, taun, Ds, h, hp, tauM, etaM, params.diffusion.mass;
                                     sources = src_mode)
    end

    # rows 1-2: consistent first moment + the correct c_M force — byte-identical
    # in structure to matrix1d_visc_HQ_BG_consistent_5f! (HQ_CM_BACKREACTION=0
    # in production keeps rows 3-5 passive)
    At2, Ax2, src2 = hqc_matrices(ϕ, tau, X[1], ur, T, dtT, drT, drur, dtur,
                                  n, dn_dT, dn_dmu, taun, Ds, h, hp; consistent = true)
    cM = get(ENV, "HQ_CM_BACKREACTION", "1") == "0" ? zero(Ds) : Ds / T
    nfloor = parse(Float64, get(ENV, "HQC_CM_NFLOOR", "1e-6"))
    cM *= n / (n + nfloor)
    at3c, ax3c, csrc = hqc_cm_row2(ϕ, tau, X[1], ur, dtur, drur, cM)

    At = zeros(MMatrix{5,5,Float64}); Ax = zeros(MMatrix{5,5,Float64})
    src = zeros(MVector{5,Float64})
    At[1, 1] = At2[1, 1]; At[1, 2] = At2[1, 2]
    At[2, 1] = At2[2, 1]; At[2, 2] = At2[2, 2]; At[2, 3] = at3c; At[2, 5] = at3c
    Ax[1, 1] = Ax2[1, 1]; Ax[1, 2] = Ax2[1, 2]
    Ax[2, 1] = Ax2[2, 1]; Ax[2, 2] = Ax2[2, 2]; Ax[2, 3] = ax3c; Ax[2, 5] = ax3c
    for j in 1:5, i in 1:3
        At[2 + i, j] = at3[i, j]
        Ax[2 + i, j] = ax3[i, j]
    end
    src[1] = src2[1]
    src[2] = src2[2] + csrc
    src[3] = src3[1]; src[4] = src3[2]; src[5] = src3[3]

    Ainv = inv(SMatrix{5,5}(At))
    A_mul_B!(A_i[1], Ainv, SMatrix{5,5}(Ax))
    jgemvavx!(Source, Ainv, SVector{5}(src))
end
