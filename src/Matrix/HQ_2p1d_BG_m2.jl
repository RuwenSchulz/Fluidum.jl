# ==============================================================================
# HQ_2p1d_BG_m2.jl — THE CONSISTENT SECOND MOMENT ON A PRESCRIBED 2+1D BACKGROUND
# (2026-09-09).
#
# The 2-D twin of HQ_const_BG_consistent_m2.jl, and the second-moment companion
# of HQ_2p1d_BG.jl (which carries the first moment only, and in its SHIPPED
# ∇α-only closure).  Fields
#
#     phi = (alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy}, piQ^eta_eta, Pi_Q)
#
# on Milne (tau, x, y, eta_s), boost invariant, metric diag(-1, 1, 1, tau^2), in
# the source-on-the-LHS convention
#
#     A_t d_tau phi + A_x d_x phi + A_y d_y phi + src = 0
#
# which is the convention of every matrix in this directory.
#
# NOT INCLUDED BY src/Fluidum.jl.  Drivers load it at runtime after
# HQ_2p1d_BG.jl and register their own mode, the same pattern
# HQ_const_BG_consistent_m2.jl uses in 1-D.
#
# ------------------------------------------------------------------------------
# WHY THE FIELDS LOOK LIKE THIS
# ------------------------------------------------------------------------------
# Boost invariance kills every eta-mixed component of piQ (odd under eta -> -eta),
# leaving tau tau, tau x, tau y, xx, xy, yy, eta eta = seven, cut by three
# orthogonality relations (u_mu piQ^{mu nu} = 0) and one trace condition to THREE
# independent dofs.  They are stored as FOUR, with piQ^eta_eta redundant and
# restored by projection — the same storage the medium shear uses in
# Matrix/2d_viscous.jl and in FiVo's src2d/shear2d.jl.  The tau-row is DERIVED
# from the transverse block, never stored, so orthogonality holds by construction.
#
# ------------------------------------------------------------------------------
# ⚠ THE CONTRACTIONS NEED THE METRIC.  THIS IS WHERE THE ERRORS LIVE.
# ------------------------------------------------------------------------------
# In the 1-D cylindrical chart the second-moment fields are components on a
# PARALLEL-TRANSPORTED ORTHONORMAL TRIAD, so every contraction looks like a plain
# product and the metric is invisible.  Here the stored fields are CONTRAVARIANT,
# and lowering an index on a u-orthogonal tensor brings in its tau-row.  Six
# distinct terms are affected, and getting any of them wrong is a silent O(u^2)
# error that looks perfectly smooth:
#
#     pi:sigma          both indices lowered  (NOT pxx sxx + 2 pxy sxy + ...)
#     pi^{(i}_l s^{j)l} one index lowered
#     a.nu, nu.grad(Th) carry -(a.u)(nu.u)/(u^tau)^2
#     the rank-1 class-(iv) block: Delta-trace over THREE directions
#     the pi:sigma trace subtraction: weighted by Delta^{ij}, not bare
#     sigma_(nu)        NOT the u-shear with u -> nu; the projector drags in
#                       gradients of U through nu^tau = (u.nu)/u^tau
#
# The FiVo twin (src2d/hq_consistent_m2_2d.jl) made all six mistakes in turn and
# its axisymmetric-reduction gate caught each one, taking it from 65 % to
# 2.3e-15.  This file is written against that record.
#
# ------------------------------------------------------------------------------
# TWO TERMS ADDED 2026-09-11 -- AND WHY NO GATE HAD CAUGHT THEM
# ------------------------------------------------------------------------------
# Both belong to the derived equation and both were simply absent, here and in
# FiVo, because both 2-D files were written from the 1-D reduction:
#
#   hq2d_proj_corr   the Delta-PROJECTION of the comoving derivative.  DEFAULT ON
#                    (`projector = true`) -- it is not optional physics, the
#                    switch exists only to reproduce a pre-2026-09-11 number.
#   hq2d_vorticity   the DNMR corotation term 2 pi^{l<i} omega_l^{j>}.  DEFAULT
#                    OFF (`vorticity = false`), to keep parity with FiVo's
#                    default and leave every existing number unchanged.
#
# ⛔ THE LESSON, and it is the same one twice: a cross-code gate proves the two
# codes agree, NOT that either is right.  N3 held at 2.2e-16 while BOTH codes
# omitted BOTH terms.  And the 1-D-reduction gates are blind by construction --
# the vorticity coupling VANISHES IDENTICALLY on an axisymmetric state, and N1
# mapped the 2-D rate onto the 1-D p_l with the projection HELD FIXED, which is
# exactly the projector term, so a frozen map agreed with projector-less rows.
# N1's mapping is corrected (it now differentiates the projection along the
# trajectory) and a new gate covers the terms on genuinely 2-D states:
# Projects/FiVoFluidumComparison/gate_2p1d_m2_newterms.jl.
#
# Both are derived symbolically, and the derivation is gated: G6 (the projector's
# closed form, exactly) and G7 (the coupling's antisymmetry, u-orthogonality,
# tracelessness, its vanishing in 1-D, and the split that fixes its SIGN --
# sigma-coupling + omega-coupling == the full pi.NablaPerp u contraction) in
# Julia/tools/derive_hq_m2_2p1d.wls, 22/22.
#
# ------------------------------------------------------------------------------
# VALIDATION
# ------------------------------------------------------------------------------
# This is an INDEPENDENT implementation, not a call into FiVo: the point of a
# cross-check is two codes, and a shared implementation would prove nothing.
#     N1  the axisymmetric reduction against Fluidum's OWN 1-D hqc_m2_rows
#     N2  rotational covariance of the traceless block
#     N3  agreement with FiVo's src2d/hq_consistent_m2_2d.jl at identical states
#     N4  the live 8x8 driver assembles those rows
#     N5  the projector moves the rows; the coupling is traceless, dies in 1-D,
#         and enters exactly as -W/u^tau            [gate_2p1d_m2.jl, 5/5]
#     M1-M5  both terms, both codes, on GENUINELY 2-D states: every switch
#         combination agrees at ~1e-16, the coupling is live there (1.4e-4), and
#         the two switches are additive  [gate_2p1d_m2_newterms.jl, 3/3]
# ==============================================================================

"""
    hq2d_m2_coeffs(T, n, taun, Ds, h, mq) -> NamedTuple

The consistent second moment's scalar coefficients — identical closed forms to
the 1-D `hqc_m2_rows` block, because they are built from the EoS and nothing
about them is chart-dependent.
"""
function hq2d_m2_coeffs(T, n, taun, Ds, h, mq)
    z  = mq / T
    K2 = besselk(2, z); K3 = besselk(3, z); K4 = besselk(4, z)
    K5 = K3 + 8 * K4 / z
    dlnC_dlnT = z * ((K3 + K5) / (2 * K4) - (K2 + K4) / (2 * K3))
    return (dlnC_dlnT = dlnC_dlnT,
            ηbar = taun * n * T / 2,              # tau_n P0/2, P0 = nT
            λa   = Ds * (mq^2 + 6 * T * h) / (2 * T),
            Aco  = mq^2 + 5 * T * h,              # A = I41/P0
            ABr  = 5 + z * K2 / K3)               # A/B = dln I31/dln T
end

"""
    hq2d_piQ_eta(ux, uy, uT, pxx, pxy, pyy)

`piQ^eta_eta` from the transverse block by tracelessness — the redundant dof.
Equivalent to `Pi^{tt} - pxx - pyy` with the tau-row derived from orthogonality.
"""
@inline function hq2d_piQ_eta(ux, uy, uT, pxx, pxy, pyy)
    ptx = (ux*pxx + uy*pxy)/uT
    pty = (ux*pxy + uy*pyy)/uT
    ptt = (ux*ptx + uy*pty)/uT
    return ptt - pxx - pyy
end

"""
    hq2d_lower2(pxx, pxy, pyy, pet, sxx, sxy, syy, set, ux, uy)

`pi_{mu nu} sigma^{mu nu}` for two transverse traceless u-orthogonal tensors
stored by their (xx, xy, yy, eta) components with the tau-row derived.

⚠ NOT `pxx sxx + 2 pxy sxy + pyy syy + pet set`.
"""
@inline function hq2d_lower2(pxx, pxy, pyy, pet, sxx, sxy, syy, set, ux, uy)
    d = 1 + ux^2 + uy^2
    return (2*pxy*sxy + pyy*syy + 2*pxy*sxy*ux^2 + 2*pyy*syy*ux^2 + pyy*syy*ux^4 -
            2*pxy*sxx*ux*uy - 2*pyy*sxy*ux*uy - 2*pxy*syy*ux*uy -
            2*pyy*sxy*ux^3*uy - 2*pxy*syy*ux^3*uy + 2*pxy*sxy*uy^2 +
            pyy*sxx*ux^2*uy^2 + 4*pxy*sxy*ux^2*uy^2 - 2*pxy*sxx*ux*uy^3 +
            pxx*sxx*(1 + uy^2)^2 + pet*set*d^2 +
            pxx*ux*uy*(syy*ux*uy - 2*sxy*(1 + uy^2))) / d^2
end

"""
    hq2d_mix(pxx, pxy, pyy, sxx, sxy, syy, ux, uy) -> (cxx, cxy, cyy)

`pi^{(i}_lambda sigma^{j)lambda}`, symmetrised, on the transverse block. One
index lowered by the metric; the tau-row contributes `-(pi.u)^i (sigma.u)^j/(u^tau)^2`.
"""
@inline function hq2d_mix(pxx, pxy, pyy, sxx, sxy, syy, ux, uy)
    d = 1 + ux^2 + uy^2
    cxx = pxx*sxx + pxy*sxy - ((pxx*ux + pxy*uy)*(sxx*ux + sxy*uy))/d
    cyy = pxy*sxy + pyy*syy - ((pxy*ux + pyy*uy)*(sxy*ux + syy*uy))/d
    cxy = (pyy*(sxy + sxy*ux^2 - sxx*ux*uy) +
           pxy*(sxx + syy + syy*ux^2 - 2*sxy*ux*uy + sxx*uy^2) +
           pxx*(sxy - syy*ux*uy + sxy*uy^2)) / (2*d)
    return (cxx, cxy, cyy)
end

"""
    hq2d_proj_corr(pxx, pxy, pyy, ax, ay, ux, uy, uT) -> (dxx, dxy, dyy)

The correction that turns a plain comoving derivative into a Delta-PROJECTED
one. The equation's first term is `tau_M Delta^{ij}_{ab} D pi_Q^{ab}`, and for a
symmetric, traceless, u-orthogonal `pi` the projection is closed-form:
differentiating `u_a pi^{ab} = 0` gives `u_a D pi^{ab} = -pi^{ab} a_a`, hence

    Delta^{ij}_{ab} D pi^{ab} = D pi^{ij} - u^i c^j - u^j c^i,   c^j = pi^{j b} a_b

so what must be SUBTRACTED from the LHS source is `tau_M (u^i c^j + u^j c^i)`.
Returned here as that quantity, already symmetrised; the caller subtracts it.

⚠ ADDED 2026-09-11. Until then this file — and FiVo's twin until 2026-09-10 —
integrated the plain `D pi^{ij}`, i.e. the correction was simply missing. The
cross-code gate could not see it, because BOTH codes omitted it and agreed at
2.2e-16 (gate_2p1d_m2.jl N3). Nor could the 1-D-reduction gate N1: it mapped the
2-D rate onto the 1-D `p_l` with the projection held FIXED (`l_mu l_nu d_tau
pi^{mu nu}`), but `l` depends on `u(tau)`, and `d_tau(l l pi) - l l d_tau pi` is
exactly this term — so the mis-mapped gate agreed with the projector-less code.

Derived and gated symbolically in `Julia/tools/derive_hq_m2_2p1d.wls` (G6):
the orthogonality identity, the closed form EXACTLY on the transverse block, and
the fact that it leaves no trace remainder. Size: 1.6e-4 to 3.4e-3 of the rate.
It vanishes at rest on the EVOLVED ij block (it is O(u a)); its tau-i components
do not, which is why only the ij block is stored.
"""
@inline function hq2d_proj_corr(pxx, pxy, pyy, ax, ay, ux, uy, uT)
    # the tau-row of pi, from orthogonality (the same construction as everywhere here)
    ptx = (ux*pxx + uy*pxy)/uT
    pty = (ux*pxy + uy*pyy)/uT
    aT  = (ux*ax + uy*ay)/uT                     # a^tau, from u.a = 0
    # c^j = pi^{j b} a_b, with a_b = g_{bb} a^b: the tau-row enters with -1
    cx = -ptx*aT + pxx*ax + pxy*ay
    cy = -pty*aT + pxy*ax + pyy*ay
    return (2*ux*cx, ux*cy + uy*cx, 2*uy*cy)
end

"""
    hq2d_vorticity(pxx, pxy, pyy, ax, ay, ux, uy, uT,
                   dtux, dtuy, dxux, dxuy, dyux, dyuy) -> (wxx, wxy, wyy)

The transverse block of the DNMR corotation term `2 pi^{l<i} omega_l^{j>}`, with

    omega^{a n} = (1/2)(NablaPerp^a u^n - NablaPerp^n u^a),
    NablaPerp^a u^n = g^{aa} d_a u^n + u^a a^n        (derivative index FIRST)

Multiply by `tau_M` and ADD to the LHS source.

⚠ `omega` IS NOT `d_x u^y - d_y u^x`. A boost-invariant transverse flow carries
the acceleration terms too: `omega^{xy} = (1/2)(d_xu^y - d_yu^x + u^x a^y - u^y a^x)`.

Every eta-component of omega vanishes by boost invariance, so only the (tau,x,y)
block enters, and

    2 pi^{l<i} omega_l^{j>} = (pi g omega) + (pi g omega)^T

is ALREADY symmetric, traceless and u-orthogonal (pi symmetric, omega
antisymmetric, both u-orthogonal) — no projection is needed. All four of those
properties are gated symbolically in `Julia/tools/derive_hq_m2_2p1d.wls` (G7).

⚠ OFF BY DEFAULT (`vorticity = false`), to preserve parity with FiVo's default
and so that every existing number is unchanged. It VANISHES IDENTICALLY in 1+1D
— a radial flow has no transverse vorticity — which is why the 1-D code never
needed it, why both 2-D ports inherited its absence from the 1-D reduction, and
why no 1-D-limit gate can see it. THE GATE THAT FIXES ITS SIGN is G7-splits:
sigma-coupling + omega-coupling equals the full contraction
`pi^{l<i} NablaPerp_l u^{j>}` EXACTLY, which ties this term's normalisation to
the sigma coupling that gate N1 already holds at 1e-15.

Measured in FiVo on a real Pb+Pb event, converged over N = 160/480/640:
`|omega|/|sigma|` median 3.01e-2 -> 2.84e-2 -> 2.79e-2, moving pi_Q by ~1e-2.
"""
@inline function hq2d_vorticity(pxx, pxy, pyy, ax, ay, ux, uy, uT,
                                dtux, dtuy, dxux, dxuy, dyux, dyuy)
    aT = (ux*ax + uy*ay)/uT                      # a^tau, from u.a = 0
    # d_i u^tau follows from u^tau = sqrt(1 + u_perp^2)
    dxuT = (ux*dxux + uy*dxuy)/uT
    dyuT = (ux*dyux + uy*dyuy)/uT
    # A^{a n} = NablaPerp^a u^n on the mixed pairs; g^{tau tau} = -1 flips the d_tau row
    Atx = -dtux + uT*ax;   Axt = dxuT + ux*aT
    Aty = -dtuy + uT*ay;   Ayt = dyuT + uy*aT
    Axy =  dxuy + ux*ay;   Ayx = dyux + uy*ax
    wtx = (Atx - Axt)/2
    wty = (Aty - Ayt)/2
    wxy = (Axy - Ayx)/2
    # pi^{tau i} from orthogonality; pi^{tau tau} does not enter the ij block
    ptx = (ux*pxx + uy*pxy)/uT
    pty = (ux*pxy + uy*pyy)/uT
    # (pi g omega)^{i n} = -pi^{i tau} omega^{tau n} + pi^{i x} omega^{x n} + pi^{i y} omega^{y n}
    Mxx = -ptx*wtx              - pxy*wxy
    Mxy = -ptx*wty + pxx*wxy
    Myx = -pty*wtx              - pyy*wxy
    Myy = -pty*wty + pxy*wxy
    return (2*Mxx, Mxy + Myx, 2*Myy)
end

"""
    hq2d_rank1(Ax, Ay, Bx, By, ux, uy) -> (qxx, qxy, qyy)

`A^{<i} B^{j>}` — symmetric, traceless, transverse — for two u-orthogonal
transverse vectors. The Delta-trace runs over THREE directions and picks up the
tau-rows of both vectors.
"""
@inline function hq2d_rank1(Ax, Ay, Bx, By, ux, uy)
    d  = 1 + ux^2 + uy^2
    tr = (Ay*(By + By*ux^2 - Bx*ux*uy) + Ax*(Bx - By*ux*uy + Bx*uy^2)) / d
    return (Ax*Bx - (1 + ux^2)*tr/3,
            (Ay*Bx + Ax*By)/2 - (ux*uy)*tr/3,
            Ay*By - (1 + uy^2)*tr/3)
end

"""
    hq2d_sigma(ux, uy, uT, tau, ax, ay, theta, dxux, dxuy, dyux, dyuy)
        -> (sxx, sxy, syy, seta)

The medium's shear `sigma^{ij}` and `sigma^eta_eta = tau^2 sigma^{eta eta}`.
"""
@inline function hq2d_sigma(ux, uy, uT, tau, ax, ay, theta, dxux, dxuy, dyux, dyuy)
    third = theta/3
    return (dxux + ux*ax - third*(1 + ux*ux),
            (dxuy + dyux)/2 + (ux*ay + uy*ax)/2 - third*(ux*uy),
            dyuy + uy*ay - third*(1 + uy*uy),
            uT/tau - third)
end

"""
    hq2d_signu(ux, uy, uT, tau, nux, nuy,
               dtnx, dxnx, dynx, dtny, dxny, dyny,
               dtux, dxux, dyux, dtuy, dxuy, dyuy) -> (sxx, sxy, syy, seta)

`sigma_(nu)^{ij}`, the transverse-projected traceless symmetric gradient of the
current.

⚠ NOT the u-shear with u -> nu. The projector acts on the FULL `grad^{(a}nu^{b)}`,
whose tau-row involves `nu^tau = (u.nu)/u^tau` — a function of u — so projecting
drags in terms carrying gradients of U. Derived symbolically
(Julia/tools/derive_signu_2p1d.wls, 3/3) and shared with the FiVo twin's
generated form, which reproduces the 1-D `(grad nu)_ll` to 1e-15.
"""
@inline function hq2d_signu(ux, uy, uT, tau, nux, nuy,
                            dtnx, dxnx, dynx, dtny, dxny, dyny,
                            dtux, dxux, dyux, dtuy, dxuy, dyuy)
    sxx = (tau*(1 + ux^2 + uy^2)*(2*dtnx*ux^3 - dtny*uy + 2*dxnx*uT - dyny*uT - ux^2*(dtny*uy - 2*dxnx*uT + dyny*uT) + ux*(2*dtnx + 3*dtnx*uy^2 + 3*dynx*uy*uT)) + nuy*((-1 - dtux*tau*ux - ux^2 + 2*dtux*tau*ux^3)*uy^3 - 2*tau*ux*(1 + ux^2)^2*(dtuy*ux + dxuy*uT) - tau*(-1 + 2*ux^2)*uy^2*(dtuy + dtuy*ux^2 - dyux*ux*uT) + (1 + ux^2)*uy*(-1 + 2*dtux*tau*ux^3 + dyuy*tau*uT + ux^2*(-1 + 2*dxux*tau*uT - 2*dyuy*tau*uT))) + nux*(ux^5*(-1 + 2*dtuy*tau*uy) + tau*uy*(1 + uy^2)*(dtux*uy + dyux*uT) + 2*tau*ux^4*(-dtux - dtux*uy^2 + dxuy*uy*uT) - 2*tau*ux^2*(dtux + 2*dtux*uy^2 + dtux*uy^4 + (-dxuy + dyux)*uy*uT + dyux*uy^3*uT) - ux*(1 + dtuy*tau*uy^3 + 2*dxux*tau*uT + uy^2*(1 + 2*dxux*tau*uT + dyuy*tau*uT)) + ux^3*(2*dtuy*tau*uy + 2*dtuy*tau*uy^3 - 2*(1 + dxux*tau*uT) + uy^2*(-1 - 2*dxux*tau*uT + 2*dyuy*tau*uT))))/(3*tau*(1 + ux^2 + uy^2)^(3/2))
    sxy = (tau*(1 + ux^2 + uy^2)*(3*dtny*ux^3 + 3*dynx*(1 + uy^2)*uT + 3*(dtnx*uy*(1 + uy^2) + dxny*uT) + ux^2*(dtnx*uy + 3*dxny*uT) + ux*(3*dtny + dtny*uy^2 + (dxnx + dyny)*uy*uT)) - nuy*(4*dtuy*tau*ux^5*uy + 3*dxuy*tau*uy*uT + 4*tau*ux^4*uy*(-(dtux*uy) + dxuy*uT) - tau*ux^2*uy*(-7*dxuy*uT + 2*dtux*uy*(3 + 2*uy^2) + dyux*uT*(3 + 4*uy^2)) + ux^3*(10*dtuy*tau*uy + 4*dtuy*tau*uy^3 + 3*dyuy*tau*uT + uy^2*(2 - 4*dxux*tau*uT + 4*dyuy*tau*uT)) + ux*(6*dtuy*tau*uy + 4*dtuy*tau*uy^3 + 2*uy^4 + 3*dyuy*tau*uT + uy^2*(2 - 3*dxux*tau*uT + 4*dyuy*tau*uT))) + nux*(2*ux^4*uy*(-1 + 2*dtuy*tau*uy) - 3*dxux*tau*uy*(1 + uy^2)*uT - 4*tau*ux^3*uy*(dtux + dtux*uy^2 - dxuy*uy*uT) + ux^2*uy*(-2 + 6*dtuy*tau*uy + 4*dtuy*tau*uy^3 - 4*dxux*tau*uT + 3*dyuy*tau*uT + uy^2*(-2 - 4*dxux*tau*uT + 4*dyuy*tau*uT)) - tau*ux*(dyux*uT*(3 + 7*uy^2 + 4*uy^4) + uy*(6*dtux + 10*dtux*uy^2 + 4*dtux*uy^4 - 3*dxuy*uy*uT))))/(6*tau*(1 + ux^2 + uy^2)^(3/2))
    syy = (tau*(1 + ux^2 + uy^2)*(-(dtnx*ux) + 2*dtny*uy + 3*dtny*ux^2*uy - dtnx*ux*uy^2 + 2*dtny*uy^3 + 3*dxny*ux*uy*uT - dxnx*(1 + uy^2)*uT + 2*dyny*(1 + uy^2)*uT) + nux*(ux^3*(-1 - dtuy*tau*uy - uy^2 + 2*dtuy*tau*uy^3) - 2*tau*uy*(1 + uy^2)^2*(dtux*uy + dyux*uT) - tau*ux^2*(-1 + 2*uy^2)*(dtux + dtux*uy^2 - dxuy*uy*uT) + ux*(1 + uy^2)*(-1 + 2*dtuy*tau*uy^3 + dxux*tau*uT + uy^2*(-1 - 2*dxux*tau*uT + 2*dyuy*tau*uT))) + nuy*((-1 + 2*dtux*tau*ux)*uy^5 + tau*ux*(1 + ux^2)*(dtuy*ux + dxuy*uT) + 2*tau*uy^4*(-dtuy - dtuy*ux^2 + dyux*ux*uT) - 2*tau*uy^2*(dtuy + 2*dtuy*ux^2 + dtuy*ux^4 + (dxuy - dyux)*ux*uT + dxuy*ux^3*uT) + uy^3*(2*dtux*tau*ux + 2*dtux*tau*ux^3 + ux^2*(-1 + 2*dxux*tau*uT - 2*dyuy*tau*uT) - 2*(1 + dyuy*tau*uT)) - uy*(1 + dtux*tau*ux^3 + 2*dyuy*tau*uT + ux^2*(1 + dxux*tau*uT + 2*dyuy*tau*uT))))/(3*tau*(1 + ux^2 + uy^2)^(3/2))
    set = (-dxnx - dyny + ((dtux*ux + dtuy*uy)*(nux*ux + nuy*uy))/(1 + ux^2 + uy^2)^(3/2) - (dtux*nux + dtuy*nuy + dtnx*ux + dtny*uy)/uT + (2*(nux*ux + nuy*uy))/(tau*uT) + nux*(dxux*ux + dyux*uy + dtux*uT) + nuy*(dxuy*ux + dyuy*uy + dtuy*uT) - ((nux*ux + nuy*uy)*(dxux*ux^2 + ux*((dxuy + dyux)*uy + dtux*uT) + uy*(dyuy*uy + dtuy*uT)))/(1 + ux^2 + uy^2))/3
    return (sxx, sxy, syy, set)
end

"""
    hq2d_m2_rows(X, tau, T, ux, uy, dtT, dxT, dyT,
                 dtux, dxux, dyux, dtuy, dxuy, dyuy,
                 dta, dtnx, dtny, n, taun, Ds, h, hp, tauM, etaM, mq)
        -> (dpxx, dpxy, dpyy, dPiQ)

`d_tau` of the transverse second-moment block and the trace, under the CONSISTENT
closure. Returned as explicit rates rather than matrix rows because rows 4-8
relax diagonally (`tau_M u^tau d_tau X + src = 0`) and carry no derivative
couplings among themselves — the same structure the 1-D `hqc_m2_rows` has, where
`at3[i, 2+i] = tau_M u^tau` and the moment block is diagonal.

`piQ^eta_eta` is NOT returned: it is redundant and is restored from the
transverse block by `hq2d_piQ_eta` after the update.

`X = (alpha, nu^x, nu^y, pxx, pxy, pyy, peta, PiQ)`.

TWO TERMS ARE SWITCHABLE, both added 2026-09-11 (see their helpers for the full
story and for why no existing gate could see either):

  * `projector = true` (DEFAULT ON) — the Delta-projection of the comoving
    derivative, `hq2d_proj_corr`. It belongs in the equation; it was simply
    missing. Turning it OFF reproduces this file's pre-2026-09-11 behaviour, and
    is the way to reproduce an old number.
  * `vorticity = false` (DEFAULT OFF) — the DNMR corotation term,
    `hq2d_vorticity`. It also belongs in the equation, but it is left off so that
    parity with FiVo's default (`m2_vorticity = false`) and every existing number
    are unchanged. Turn BOTH codes on together.
"""
function hq2d_m2_rows(X, tau, T, ux, uy, dtT, dxT, dyT,
                      dtux, dxux, dyux, dtuy, dxuy, dyuy,
                      dxa, dya, dta, dtnx, dxnx, dynx, dtny, dxny, dyny,
                      n, taun, Ds, h, hp, tauM, etaM, mq;
                      projector::Bool = true, vorticity::Bool = false)
    uT = sqrt(1 + ux^2 + uy^2)
    nux = X[2]; nuy = X[3]
    pxx = X[4]; pxy = X[5]; pyy = X[6]; pet = X[7]; PiQ = X[8]
    tauM > 0 || return (zero(T), zero(T), zero(T), zero(T))

    ηQ = etaM
    ζQ = (5//3) * ηQ
    cf = hq2d_m2_coeffs(T, n, taun, Ds, h, mq)

    # kinematics of the background
    ax = uT*dtux + ux*dxux + uy*dyux
    ay = uT*dtuy + ux*dxuy + uy*dyuy
    dtuT  = (ux*dtux + uy*dtuy)/uT
    theta = dtuT + dxux + dyuy + uT/tau

    sxx, sxy, syy, set = hq2d_sigma(ux, uy, uT, tau, ax, ay, theta,
                                    dxux, dxuy, dyux, dyuy)
    nxx, nxy, nyy, _   = hq2d_signu(ux, uy, uT, tau, nux, nuy,
                                    dtnx, dxnx, dynx, dtny, dxny, dyny,
                                    dtux, dxux, dyux, dtuy, dxuy, dyuy)

    d    = 1 + ux^2 + uy^2
    DT   = uT*dtT + ux*dxT + uy*dyT
    DlnT = DT / T
    DlnC = cf.dlnC_dlnT * DlnT
    Dα   = uT*dta + ux*dxa + uy*dya
    geo  = tauM * ((5//3)*theta + DlnC)

    πσ   = hq2d_lower2(pxx, pxy, pyy, pet, sxx, sxy, syy, set, ux, uy)
    cxx, cxy, cyy = hq2d_mix(pxx, pxy, pyy, sxx, sxy, syy, ux, uy)

    # theta_nu and a.nu, both with the metric
    nt   = (ux*nux + uy*nuy)/uT
    dtnt = (ux*dtnx + uy*dtny + nux*dtux + nuy*dtuy)/uT - nt*(ux*dtux + uy*dtuy)/(uT*uT)
    thnu = dtnt + dxnx + dyny + nt/tau
    adn  = ax*nux + ay*nuy - ((ax*ux + ay*uy)*(nux*ux + nuy*uy))/d

    # class (iv)
    gTx = (h + T*hp) * (dxT + ux*DT)
    gTy = (h + T*hp) * (dyT + uy*DT)
    qa  = hq2d_rank1(ax, ay, nux, nuy, ux, uy)
    qg  = hq2d_rank1(nux, nuy, gTx, gTy, ux, uy)
    qxx = 2*cf.λa*qa[1] + (Ds/T)*qg[1]
    qxy = 2*cf.λa*qa[2] + (Ds/T)*qg[2]
    qyy = 2*cf.λa*qa[3] + (Ds/T)*qg[3]

    # the traceless LHS sources
    Sxx = pxx + 2*ηQ*nxx + geo*pxx + 2*tauM*(cxx - (1 + ux*ux)*πσ/3) +
          2*tauM*PiQ*sxx + 2*cf.ηbar*sxx + qxx
    Syy = pyy + 2*ηQ*nyy + geo*pyy + 2*tauM*(cyy - (1 + uy*uy)*πσ/3) +
          2*tauM*PiQ*syy + 2*cf.ηbar*syy + qyy
    Sxy = pxy + 2*ηQ*nxy + geo*pxy + 2*tauM*(cxy - (ux*uy)*πσ/3) +
          2*tauM*PiQ*sxy + 2*cf.ηbar*sxy + qxy

    # the corotation (vorticity) coupling — ADDED to the LHS source, see hq2d_vorticity
    if vorticity
        wxx, wxy, wyy = hq2d_vorticity(pxx, pxy, pyy, ax, ay, ux, uy, uT,
                                       dtux, dtuy, dxux, dxuy, dyux, dyuy)
        Sxx += tauM*wxx; Sxy += tauM*wxy; Syy += tauM*wyy
    end

    # the Delta-projection of D pi — SUBTRACTED, see hq2d_proj_corr
    if projector
        dxx, dxy, dyy = hq2d_proj_corr(pxx, pxy, pyy, ax, ay, ux, uy, uT)
        Sxx -= tauM*dxx; Sxy -= tauM*dxy; Syy -= tauM*dyy
    end

    # the trace channel — a scalar, so identical in form to 1-D
    ngT = nux*gTx + nuy*gTy - ((nux*ux + nuy*uy)*(gTx*ux + gTy*uy))/d
    SB  = PiQ + ζQ*thnu + geo*PiQ + (2//3)*tauM*πσ +
          cf.ηbar*(cf.ABr*DlnT + (5//3)*theta) +
          (Ds*cf.Aco/(3*T))*adn + (5//6)*(Ds/T)*ngT + cf.ηbar*Dα

    den = tauM * uT
    return (-Sxx/den, -Sxy/den, -Syy/den, -SB/den)
end

# ==============================================================================
# THE LIVE DRIVER
# ==============================================================================

"""
    matrix2d_HQ_BG_m2!(; dmn_eps = 1e-6, background_fields, α_max = 200.0, consistent_fm = true)

Closure constructor for the 8-field 2+1D charm sector WITH the consistent second
moment, matching `matrix2d_HQ_BG!`'s signature. Returns the
`(A_i, Source, ϕ, tau, X, params)` callable `basicupwinding` expects, with
`X = (x, y)` and `A_i` a 2-tuple of 8×8 matrices.

`background_fields(tau, x, y)` returns the same 12-tuple `matrix2d_HQ_BG!` takes:
`(T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy)`.

⚠ ROWS 4-8 ARE PASSIVE AT c_M = 0, WHICH IS THE ONLY CASE THIS DRIVER SUPPORTS.
The second moments read the current and the background but feed back into
neither, so the 8×8 is BLOCK TRIANGULAR: rows 1-3 are exactly
`matrix2d_HQ_BG!`'s, unchanged, and rows 4-8 hang off them. That is what makes it
safe to reuse the first-moment block verbatim rather than re-deriving it, and it
is the same structure the 1-D `matrix1d_visc_HQ_BG_consistent_m2_5f!` has with
`HQ_CM_BACKREACTION=0`.
"""
function matrix2d_HQ_BG_m2!(; dmn_eps = 1e-6, background_fields, α_max = 200.0,
                              consistent_fm::Bool = true,
                              projector::Bool = true, vorticity::Bool = false)
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_HQ_BG_m2!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields,
        α_max = α_max, consistent_fm = consistent_fm,
        projector = projector, vorticity = vorticity)
end

function matrix2d_HQ_BG_m2!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0,
        consistent_fm::Bool = true,
        projector::Bool = true, vorticity::Bool = false)

    T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy =
        background_fields(tau, X[1], X[2])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    uT = sqrt(1 + ux^2 + uy^2)

    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_da = thermo.pressure_derivative
    dn_da += dmn_eps

    κ    = diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    taun = τ_diffusion_hadron(T, α_safe, params.eos, params.diffusion) * HQ_TAUN_SCALE[]

    # ── rows 1-3: the first moment, VERBATIM from the 3-field driver ──────────
    # h, h' and D_s are needed by BOTH the consistent first moment and the second
    # moment, so they are built before either.
    h, taun_h, Ds = hqc_enthalpy_per_particle(T, α_safe, params.eos, params.diffusion)
    δ  = max(1e-4, 1e-3 * T)
    hp = (hqc_enthalpy_per_particle(T + δ, α_safe, params.eos, params.diffusion)[1] -
          hqc_enthalpy_per_particle(T - δ, α_safe, params.eos, params.diffusion)[1]) / (2δ)

    # 🔴 2026-09-09. `consistent_fm` used to be accepted and NEVER READ -- the name
    # appeared in this signature and nowhere in the body -- so every 2-D run was
    # silently on the SHIPPED (grad-alpha-only) row whatever the caller asked for.
    # Measured against a closed-form at-rest reference, the shipped row returns 0.0
    # for the temperature-gradient channel, and on a real Pb+Pb cell that channel is
    # 1.8x the alpha-drive: the DOMINANT term. See
    # Projects/FiVoFluidumComparison/COMPARISON_2P1D.md section 12.
    #
    # `hq2d_matrices_consistent` is the full-grad-P row,
    #     (D_s/T) Delta^mu_lambda Nabla_nu T_Q^{nu lambda} + nu^mu = 0,
    # generated by Julia/tools/derive_hq_consistent_2p1d.wls (15/15 gates: the
    # cylindrical reduction reproduces Fluidum's own 1-D consistent source, and the
    # at-rest limit is the closed form to 2.2e-16).
    At3, Ax3, Ay3, src3 = consistent_fm ?
        hq2d_matrices_consistent(α_safe, ϕ[2], ϕ[3], tau, T, ux, uy,
                                 dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy,
                                 n, dn_dT, dn_da, taun, κ, h, hp, Ds) :
        hq2d_matrices(α_safe, ϕ[2], ϕ[3], tau, T, ux, uy,
                      dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy,
                      n, dn_dT, dn_da, taun, κ)

    # ── rows 4-8: the consistent second moment ────────────────────────────────
    # They relax DIAGONALLY (tau_M u^tau d_tau X + tau_M u^i d_i X + src = 0) and
    # carry no derivative couplings among themselves, so the block is
    # At[4:8,4:8] = tau_M u^tau I, Ax = tau_M u^x I, Ay = tau_M u^y I, and every
    # other dependence -- on alpha, on nu, on the background -- lands in `src`.
    z    = params.diffusion.mass / T
    tauM = taun * besselk(4, z) * besselk(2, z) / (2 * besselk(3, z)^2)
    etaM = T * taun / 2

    # ⚠ THE alpha AND nu GRADIENTS GO IN THE MATRIX, NOT THE SOURCE. Rows 4-8
    # depend on d_tau alpha, d_i alpha and d nu -- in this convention those are
    # COLUMNS of At/Ax/Ay (columns 1, 2, 3), because alpha, nu^x and nu^y are
    # evolved fields of this very system. Putting them in `src` with their
    # current values frozen would lag them by a step and, worse, hide a real
    # coupling from the characteristic analysis.
    #
    # `hq2d_m2_rows` returns explicit RATES, so the columns are read off by
    # finite-differencing it in each gradient. They are exactly linear in every
    # one of them (the system is quasi-linear -- gate G3 of
    # derive_hq_m2_2p1d.wls), so a single difference is exact, not an
    # approximation.
    m2rate(dta, dxa, dya, dtnx, dxnx, dynx, dtny, dxny, dyny) =
        hq2d_m2_rows((α_safe, ϕ[2], ϕ[3], ϕ[4], ϕ[5], ϕ[6], ϕ[7], ϕ[8]),
                     tau, T, ux, uy, dtT, dxT, dyT,
                     dtux, dxux, dyux, dtuy, dxuy, dyuy,
                     dxa, dya, dta, dtnx, dxnx, dynx, dtny, dxny, dyny,
                     n, taun, Ds, h, hp, tauM, etaM, params.diffusion.mass;
                     projector = projector, vorticity = vorticity)

    Z = ntuple(_ -> zero(T), 9)
    r0 = m2rate(Z...)
    unit(k) = ntuple(i -> i == k ? one(T) : zero(T), 9)
    cols = ntuple(k -> begin
        rk = m2rate(unit(k)...)
        (rk[1] - r0[1], rk[2] - r0[2], rk[3] - r0[3], rk[4] - r0[4])
    end, 9)

    At = zeros(MMatrix{8,8,Float64}); Ax = zeros(MMatrix{8,8,Float64})
    Ay = zeros(MMatrix{8,8,Float64}); src = zeros(MVector{8,Float64})
    for i in 1:3, j in 1:3
        At[i,j] = At3[i,j]; Ax[i,j] = Ax3[i,j]; Ay[i,j] = Ay3[i,j]
    end
    src[1] = src3[1]; src[2] = src3[2]; src[3] = src3[3]
    for k in 4:8
        At[k,k] = tauM * uT; Ax[k,k] = tauM * ux; Ay[k,k] = tauM * uy
    end
    # LHS sources: At[k,k] * (explicit rate) = -src[k]  =>  src[k] = -At[k,k]*rate.
    # Row 7 (piQ^eta_eta) is the REDUNDANT dof: it is not evolved independently,
    # it is reconstructed from the transverse block by tracelessness. Relaxing it
    # toward that value keeps it consistent without adding a constraint equation.
    # rows 4-8 at zero (alpha, nu) gradients -> the source; the gradient
    # dependence -> the matrix columns. `cols[k]` is the derivative of the four
    # rates w.r.t. gradient k, ordered
    #   1: d_tau alpha  2: d_x alpha  3: d_y alpha
    #   4: d_tau nu^x   5: d_x nu^x   6: d_y nu^x
    #   7: d_tau nu^y   8: d_x nu^y   9: d_y nu^y
    ROW = (4, 5, 6, 8)                     # rows 4,5,6 traceless; 8 the trace
    for (a, k) in enumerate(ROW)
        src[k] = -(tauM*uT) * r0[a]
        At[k,1] = -(tauM*uT) * cols[1][a]; Ax[k,1] = -(tauM*uT) * cols[2][a]
        Ay[k,1] = -(tauM*uT) * cols[3][a]
        At[k,2] = -(tauM*uT) * cols[4][a]; Ax[k,2] = -(tauM*uT) * cols[5][a]
        Ay[k,2] = -(tauM*uT) * cols[6][a]
        At[k,3] = -(tauM*uT) * cols[7][a]; Ax[k,3] = -(tauM*uT) * cols[8][a]
        Ay[k,3] = -(tauM*uT) * cols[9][a]
    end
    # Row 7 (piQ^eta_eta) is the REDUNDANT dof: not evolved independently, but
    # reconstructed from the transverse block by tracelessness. Relaxing it toward
    # that value keeps it consistent without adding a constraint equation.
    # 🔴 2026-09-09. This used to read
    #     src[7] = phi[7] - hq2d_piQ_eta(...);  At[7,7] = 1
    # i.e. d_tau piQ^eta_eta = -(piQ^eta_eta - target), a relaxation with a time
    # constant of ONE fm/c. Tracelessness is an ALGEBRAIC constraint, not a slow
    # mode: over a 0.4 fm run that relaxation cannot reach its target, and measured
    # on the physical IC it did not -- piQ^eta/(pxx+pyy) came out -0.343 against the
    # -1 the constraint demands (FiVo, which imposes it, gives -0.995). The lag
    # feeds straight back through pi:sigma and the trace, and it is what froze the
    # cross-code difference at ~40% in pi_Q^xx and ~100% in Pi_Q REGARDLESS OF
    # RESOLUTION -- the signature of a constraint violation rather than a
    # discretisation error.
    #
    # `RELAX_RATE` sets how fast the constraint is enforced. It must be fast
    # compared with tau_M (the physical relaxation of the moments themselves), so
    # that row 7 tracks the transverse block instead of lagging it. 1/(0.05 fm) is
    # ~25x faster than a typical tau_M ~ 0.25 fm here and still finite, which keeps
    # the system an ODE rather than a DAE.
    RELAX_RATE = 100.0
    src[7] = RELAX_RATE * (ϕ[7] - hq2d_piQ_eta(ux, uy, uT, ϕ[4], ϕ[5], ϕ[6]))
    At[7,7] = one(T); Ax[7,7] = zero(T); Ay[7,7] = zero(T)

    Ainv = inv(SMatrix(At))
    if !all(isfinite, Ainv)
        @warn "matrix2d_HQ_BG_m2!: non-finite inv(A_t)" tau x=X[1] y=X[2] T ux uy tauM maxlog=5
    end
    A_mul_B!(A_i[1], Ainv, SMatrix(Ax))
    A_mul_B!(A_i[2], Ainv, SMatrix(Ay))
    jgemvavx!(Source, Ainv, SVector(src))
end
