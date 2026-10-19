# ==============================================================================
# Thermodynamically consistent charm second moment on a prescribed 1+1D
# cylindrical background.
#
# The background (T, u) is prescribed, not evolved. `background_fields(tau, r)`
# returns (T, ur, dtT, drT, drur, dtur), as in HQ_const_BG.jl. This file uses
# hqc_first_moment_block, hqc_matrices and hqc_cm_row2 from
# HQ_const_BG_consistent.jl.
#
# Rows 1-2 (alpha, nu^r) are the same as in matrix1d_visc_HQ_BG_consistent_5f!
# (one shared body, hqc_first_moment_block). Rows 3-5 are the thermodynamically
# consistent second-moment equations in the delta-M-extension convention (which
# agrees better with kinetic simulations than the strict-Pi variant), reduced
# for boost-invariant, azimuthally symmetric Milne flow, u = (u^tau, u^r, 0, 0),
# mostly-plus g = diag(-1, 1, r^2, tau^2).
#
# Field convention (the same as the c_M coupling): the transverse basis
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
#     tau_M D p_i + p_i + 2 eta_Q signu_i                              (i)
#       + tau_M [ ((5/3)theta + D ln C) p_i + 2(p_i s_i - (pi:sigma)/3)
#                 + 2 Pi_Q s_i ]                                       (iii + dM-ext)
#       + 2 etabar s_i                                                 (ii)
#       + [ 2 lambda_a a_l nu_l + (D_s/T) nu_l D_l(Th) ] (delta_il - 1/3)  (iv)
#     = 0,          pi:sigma = p_l s_l + p_phi s_phi + p_eta s_eta
#
# Trace channel:
#     tau_M D Pi_Q + Pi_Q + zeta_Q theta_nu                            (i)
#       + tau_M ((5/3)theta + D ln C) Pi_Q + (2/3) tau_M pi:sigma      (iii + dM-ext)
#       + (tau_n P0/2) [ D alpha + (A/B) D ln T + (5/3) theta ]        (ii)
#       + (D_s A/3T) a.nu + (5 D_s/6T) nu_l D_l(Th)                    (iv)
#     = 0
#
# Class (i) is the nabla-nu-only relaxation; (ii)-(iv) are the completion terms.
#
# COEFFICIENTS (z = m_q/T, K_n = besselk):
#     tau_M  = tau_n K4 K2 / (2 K3^2)        eta_Q  = T tau_n / 2   (= etaM)
#     zeta_Q = (5/3) eta_Q                   etabar = tau_n P0 / 2,  P0 = n T
#     lambda_a = (D_s/2T)(m^2 + 6 T h)       A = m^2 + 5 T h,  A/B = 5 + z K2/K3
#     D ln C = [z ((K3+K5)/(2K4) - (K2+K4)/(2K3))] D ln T,   K5 = K3 + 8 K4/z
#     grad(Th) = (h + T h') grad T           (h EoS-agnostic, as in rows 1-2)
# For a multi-species hadron list, h and n are the list's aggregates. A, A/B, C
# and D ln C are moments of the charm distribution, so for a list they are its
# species sums; `hq_m2_mixture` (EquationOfState/heavyquark.jl) supplies them, as
# for the two-way system. The single-quark Bessel forms at diffusion.mass differ
# from these sums by up to ~37 % for a 77-state charm hadron list. Keyword `mixture`
# (env FLUIDUM_HQ_M2_MIXTURE=0) selects the single-quark forms; a one-species
# list of mass diffusion.mass gets them in either case.
#
# The source switch: env LP1_M2_SOURCES, read once in Fluidum.__init__ into
# HQ_M2_SOURCES[], overridable per closure with the keyword `sources`:
#
#   "consistent" (default)  the equations above, all completion terms on.
#   "legacy"                rows 3-5 delegated verbatim to the older rows of
#                           one_d_viscous_matrix_fugacity_BG_only_second_moment
#                           (cM = 0), i.e. identical to
#                           matrix1d_visc_HQ_BG_consistent_5f!. Reproduces
#                           earlier results.
#   "nabla_nu"              the consistent rows with the completion terms zeroed,
#                           i.e. the nabla-nu-only relaxation
#                           tau_M D pi + pi = -2 eta_Q sigma_(nu),
#                           tau_M D Pi + Pi = -zeta_Q theta_nu
#                           in the covariant reduction.
#
# "legacy" is not the nabla_nu limit. In the older rows (1) the sigma_(nu) drive
# has the opposite sign, and (2) there are coordinate-transport couplings among
# (pi_l, pi_phi, Pi_Q) that are absent from the Delta-projected D pi in the
# parallel-transported basis. The consistent-vs-legacy difference therefore
# contains both repairs on top of classes ii-iv; nabla_nu separates them.
#
# c_M: rows 3-5 carry no cM entries of their own and are passive at the
# default cM = 0; the row-2 coupling is hqc_cm_row2, as in the first-moment file.
#
# The general-r reduction agrees pointwise with the covariant assembly, and a
# manufactured-field finite-difference residual is 1.5e-11.
# ==============================================================================

function _parse_m2_sources(s::AbstractString)
    s == "consistent" && return :consistent
    s == "legacy" && return :legacy
    s == "nabla_nu" && return :nabla_nu
    error("LP1_M2_SOURCES must be \"consistent\", \"legacy\" or \"nabla_nu\", got \"$s\"")
end
function _parse_m2_sources(s::Symbol)
    s in (:consistent, :legacy, :nabla_nu) ||
        error("LP1_M2_SOURCES / sources must be :consistent, :legacy or :nabla_nu, got :$s")
    return s
end

"""
    hqc_m2_sources() -> Symbol

The LP1_M2_SOURCES switch as this process read it (Fluidum.__init__):
:consistent (default), :legacy (the older rows verbatim) or :nabla_nu
(consistent minus the completion terms). Anything else is an error at load
time, so a misspelled value cannot silently run different physics. Assign
`Fluidum.HQ_M2_SOURCES[]` (or pass `sources =` to the closure constructor) to
change it in-process.
"""
hqc_m2_sources() = HQ_M2_SOURCES[]

"""
    hqc_m2_rows(X, tau, r, ur, T, dtT, drT, drur, dtur,
                n, taun, Ds, h, hp, tauM, etaM, mq; sources = :consistent, mix = nothing)

Rows 3-5 of the consistent-second-moment system: returns
(at3, ax3, src3) as 3x5 SMatrix rows / length-3 source, in the
source-on-the-LHS convention. `sources = :nabla_nu` zeroes the completion
terms; `:legacy` is handled by the caller (the older rows, see the file header).

`mix = nothing` (default): the EOS coefficients A/B, d ln C/d ln T and the m² of
A = I41/P0 are the single-quark closed forms at m = `mq`.
`mix = hq_m2_mixture(T, α, list, mq)` (what the driver passes for a
multi-species list): the list's species sums, and the D α coefficient becomes
η̄·∂_α ln I31 — the same coefficients as the two-way wrappers. `tauM` is the
caller's (the driver multiplies it by `mix.Cratio`).
"""
function hqc_m2_rows(X, tau, r, ur, T, dtT, drT, drur, dtur,
                     n, taun, Ds, h, hp, tauM, etaM, mq;
                     sources::Symbol = :consistent, mix = nothing)
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

    # ── class (i) rows: tau_M D p_i + p_i = -2 eta_Q signu_i ; trace with zeta_Q theta_nu
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
        if mix === nothing      # single-quark closed forms at m = mq
            z  = mq / T
            K2 = besselk(2, z); K3 = besselk(3, z); K4 = besselk(4, z)
            K5 = K3 + 8 * K4 / z
            dlnC_dlnT = z * ((K3 + K5) / (2 * K4) - (K2 + K4) / (2 * K3))
            m2  = mq^2
            ABr = 5 + z * K2 / K3                        # A/B = dln I31/dln T
            aα  = one(T)                                 # d_alpha ln I31
        else                    # the list's species sums (hq_m2_mixture)
            dlnC_dlnT = mix.dlnC_dlnT
            m2  = mix.m2bar
            ABr = mix.ABr
            aα  = mix.aI31
        end
        DlnT = (uτ * dtT + ur * drT) / T
        DlnC = dlnC_dlnT * DlnT
        DlTh = (h + T * hp) * (ur * dtT + uτ * drT)     # D_l(Th)

        ηbar = taun * n * T / 2                          # tau_n P0 / 2
        λa   = Ds * (m2 + 6 * T * h) / (2 * T)           # (Ds/2T)(A+B)
        Aco  = m2 + 5 * T * h                            # A = I41/P0
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
        if mix === nothing
            at3[3, 1] = ηbar * uτ
            ax3[3, 1] = ηbar * ur
        else
            at3[3, 1] = ηbar * aα * uτ
            ax3[3, 1] = ηbar * aα * ur
        end
    end

    return SMatrix{3,5}(at3), SMatrix{3,5}(ax3), SVector{3}(src)
end

function matrix1d_visc_HQ_BG_consistent_m2_5f!(; dmn_eps = 1e-6, background_fields,
                                               α_max = 200.0, use_NR_tauM::Bool = false,
                                               sources::Symbol = HQ_M2_SOURCES[],
                                               mixture::Bool = HQ_M2_MIXTURE[])
    src_mode = _parse_m2_sources(sources)
    (A_i, Source, ϕ, tau, X, params) -> matrix1d_visc_HQ_BG_consistent_m2_5f!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields,
        α_max = α_max, use_NR_tauM = use_NR_tauM, sources = src_mode, mixture = mixture)
end

# `mixture` (default `HQ_M2_MIXTURE[]`, env FLUIDUM_HQ_M2_MIXTURE, true): for a multi-species charm list
# the second-moment EOS coefficients (τ_M ∝ C, d ln C/d ln T, A/B, A) are the list's species sums
# (`hq_m2_mixture`, heavyquark.jl), as in the two-way system. false, a one-species list of mass
# diffusion.mass, or use_NR_tauM: the single-quark forms. Not used by :legacy.
function matrix1d_visc_HQ_BG_consistent_m2_5f!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0,
        use_NR_tauM::Bool = false, sources::Symbol = HQ_M2_SOURCES[],
        mixture::Bool = HQ_M2_MIXTURE[])

    _parse_m2_sources(sources)       # reject a misspelled mode
    # rows 1-2 and the c_M row-2 pieces: the same body as matrix1d_visc_HQ_BG_consistent_5f!
    b  = hqc_first_moment_block(ϕ, tau, X[1], params, background_fields(tau, X[1]),
                                dmn_eps, α_max, use_NR_tauM)
    tr = b.tr

    # rows 3-5: the consistent second moment; under :legacy the older rows verbatim
    # (see the header on why they are not the nabla_nu limit); under :nabla_nu the
    # covariant nabla-nu-only forms.
    local at3, ax3, src3
    if sources == :legacy
        Ats, Axs, srcs = one_d_viscous_matrix_fugacity_BG_only_second_moment(
            ϕ, tau, X[1], b.ur, b.T, b.dtT, b.drT, b.drur, b.dtur,
            b.n, b.dn_dmu, b.dn_dT, tr.taun, tr.Ds * b.n, tr.tauM, tr.etaM,
            params.diffusion.mass, zero(tr.Ds))
        at3 = SMatrix{3,5,Float64}(ntuple(k -> Ats[3 + (k - 1) % 3, 1 + (k - 1) ÷ 3], Val(15)))
        ax3 = SMatrix{3,5,Float64}(ntuple(k -> Axs[3 + (k - 1) % 3, 1 + (k - 1) ÷ 3], Val(15)))
        src3 = SVector{3,Float64}(srcs[3], srcs[4], srcs[5])
    elseif !mixture || use_NR_tauM || hq_m2_single(params.eos.hadron_list, params.diffusion.mass)
        at3, ax3, src3 = hqc_m2_rows(ϕ, tau, X[1], b.ur, b.T, b.dtT, b.drT, b.drur, b.dtur,
                                     b.n, tr.taun, tr.Ds, tr.h, b.hp, tr.tauM, tr.etaM,
                                     params.diffusion.mass; sources = sources)
    else
        mix = hq_m2_mixture(b.T, clamp(ϕ[1], -α_max, α_max), params.eos.hadron_list,
                            params.diffusion.mass)
        at3, ax3, src3 = hqc_m2_rows(ϕ, tau, X[1], b.ur, b.T, b.dtT, b.drT, b.drur, b.dtur,
                                     b.n, tr.taun, tr.Ds, tr.h, b.hp, tr.tauM * mix.Cratio, tr.etaM,
                                     params.diffusion.mass; sources = sources, mix = mix)
    end

    At2, Ax2, src2 = b.At2, b.Ax2, b.src2
    At = zeros(MMatrix{5,5,Float64}); Ax = zeros(MMatrix{5,5,Float64})
    src = zeros(MVector{5,Float64})
    At[1, 1] = At2[1, 1]; At[1, 2] = At2[1, 2]
    At[2, 1] = At2[2, 1]; At[2, 2] = At2[2, 2]; At[2, 3] = b.atc; At[2, 5] = b.atc
    Ax[1, 1] = Ax2[1, 1]; Ax[1, 2] = Ax2[1, 2]
    Ax[2, 1] = Ax2[2, 1]; Ax[2, 2] = Ax2[2, 2]; Ax[2, 3] = b.axc; Ax[2, 5] = b.axc
    for j in 1:5, i in 1:3
        At[2 + i, j] = at3[i, j]
        Ax[2 + i, j] = ax3[i, j]
    end
    src[1] = src2[1]
    src[2] = src2[2] + b.csrc
    src[3] = src3[1]; src[4] = src3[2]; src[5] = src3[3]

    Ainv = inv(SMatrix{5,5}(At))
    A_mul_B!(A_i[1], Ainv, SMatrix{5,5}(Ax))
    jgemvavx!(Source, Ainv, SVector{5}(src))
end
