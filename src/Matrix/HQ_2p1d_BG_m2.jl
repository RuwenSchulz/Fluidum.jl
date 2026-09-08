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
# VALIDATION
# ------------------------------------------------------------------------------
# This is an INDEPENDENT implementation, not a call into FiVo: the point of a
# cross-check is two codes, and a shared implementation would prove nothing.  It
# is gated three ways, in Projects/FiVoFluidumComparison/converge_2p1d_m2.jl:
#     N1  the axisymmetric reduction against Fluidum's OWN 1-D hqc_m2_rows
#     N2  rotational covariance of the traceless block
#     N3  agreement with FiVo's src2d/hq_consistent_m2_2d.jl at identical states
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
"""
function hq2d_m2_rows(X, tau, T, ux, uy, dtT, dxT, dyT,
                      dtux, dxux, dyux, dtuy, dxuy, dyuy,
                      dxa, dya, dta, dtnx, dxnx, dynx, dtny, dxny, dyny,
                      n, taun, Ds, h, hp, tauM, etaM, mq)
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

    # the trace channel — a scalar, so identical in form to 1-D
    ngT = nux*gTx + nuy*gTy - ((nux*ux + nuy*uy)*(gTx*ux + gTy*uy))/d
    SB  = PiQ + ζQ*thnu + geo*PiQ + (2//3)*tauM*πσ +
          cf.ηbar*(cf.ABr*DlnT + (5//3)*theta) +
          (Ds*cf.Aco/(3*T))*adn + (5//6)*(Ds/T)*ngT + cf.ηbar*Dα

    den = tauM * uT
    return (-Sxx/den, -Sxy/den, -Syy/den, -SB/den)
end
