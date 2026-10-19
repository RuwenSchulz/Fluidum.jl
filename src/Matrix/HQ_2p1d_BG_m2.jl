# ==============================================================================
# Consistent charm second moment on a prescribed 2+1D background.
#
# The 2-D analogue of HQ_const_BG_consistent_m2.jl, and the second-moment
# companion of HQ_2p1d_BG.jl (first moment only, grad-alpha-only closure).
# The background (T, u) is prescribed, not evolved. Fields
#
#     phi = (alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy}, piQ^eta_eta, Pi_Q)
#
# on Milne (tau, x, y, eta_s), boost invariant, metric diag(-1, 1, 1, tau^2), in
# the source-on-the-LHS convention
#
#     A_t d_tau phi + A_x d_x phi + A_y d_y phi + src = 0
#
# as everywhere in this directory. Uses `hq_transport`
# (EquationOfState/heavyquark.jl), `hqc_dh_dT` / `hqc_enthalpy_per_particle`
# (HQ_const_BG_consistent.jl) and `hq2d_matrices[_consistent]` (the two
# generated 2-D files).
#
# ------------------------------------------------------------------------------
# Field layout
# ------------------------------------------------------------------------------
# Boost invariance kills every eta-mixed component of piQ (odd under eta -> -eta),
# leaving tau tau, tau x, tau y, xx, xy, yy, eta eta = seven, cut by three
# orthogonality relations (u_mu piQ^{mu nu} = 0) and one trace condition to three
# independent dofs. They are stored as four, with piQ^eta_eta redundant and
# held on the tracelessness constraint by a stiff relaxation in the RHS (row 7 of
# the driver), the same storage the medium shear uses in Matrix/2d_viscous.jl.
# The tau-row is derived from the transverse block, never stored, so
# orthogonality holds by construction.
#
# ------------------------------------------------------------------------------
# ⚠ The contractions need the metric
# ------------------------------------------------------------------------------
# In the 1-D cylindrical chart the second-moment fields are components on a
# parallel-transported orthonormal triad, so every contraction looks like a plain
# product and the metric is invisible. Here the stored fields are contravariant,
# and lowering an index on a u-orthogonal tensor brings in its tau-row. Six
# terms are affected, and getting any of them wrong is a silent O(u^2) error
# that looks smooth:
#
#     pi:sigma          both indices lowered  (NOT pxx sxx + 2 pxy sxy + ...)
#     pi^{(i}_l s^{j)l} one index lowered
#     a.nu, nu.grad(Th) carry -(a.u)(nu.u)/(u^tau)^2
#     the rank-1 class-(iv) block: Delta-trace over THREE directions
#     the pi:sigma trace subtraction: weighted by Delta^{ij}, not bare
#     sigma_(nu)        NOT the u-shear with u -> nu; the projector drags in
#                       gradients of U through nu^tau = (u.nu)/u^tau
#
# ------------------------------------------------------------------------------
# Two terms absent from the 1-D reduction
# ------------------------------------------------------------------------------
#   hq2d_proj_corr   the Delta-projection of the comoving derivative. On by
#                    default (`projector = true`); `false` reproduces earlier
#                    results.
#   hq2d_vorticity   the DNMR corotation term 2 pi^{l<i} omega_l^{j>}. Off by
#                    default (`vorticity = false`), so existing results are
#                    unchanged.
#
# Neither shows up in a 1-D reduction: the vorticity coupling vanishes
# identically on an axisymmetric state, and the projector term is lost if the
# 2-D rate is mapped onto the 1-D p_l with the projection held fixed. Both are
# derived symbolically: the projector's closed form, and the coupling's
# antisymmetry, u-orthogonality, tracelessness, vanishing in 1-D and sign
# (sigma-coupling + omega-coupling = the full pi.NablaPerp u contraction).
#
# ------------------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------------------
#     the axisymmetric reduction reproduces the 1-D hqc_m2_rows;
#     the traceless block is rotationally covariant;
#     an independent finite-volume implementation agrees at identical states to
#         ~1e-16, including genuinely 2-D states with every switch combination;
#     the 8x8 driver assembles those rows.
# ==============================================================================

"""
    hq2d_m2_coeffs(T, n, taun, Ds, h, mq) -> NamedTuple

The consistent second moment's scalar coefficients, the same closed forms as in
the 1-D `hqc_m2_rows` (they come from the EoS and are chart-independent). These
are the single-quark forms at m = `mq`; for a multi-species list the driver uses
the method below. For a 77-state charm hadron list the single-quark A/B,
d ln C/d ln T, C and A differ from the species sums by 2–37 %.
"""
function hq2d_m2_coeffs(T, n, taun, Ds, h, mq)
    z  = mq / T
    K2 = besselk(2, z); K3 = besselk(3, z); K4 = besselk(4, z)
    K5 = K3 + 8 * K4 / z
    dlnC_dlnT = z * ((K3 + K5) / (2 * K4) - (K2 + K4) / (2 * K3))
    ηbar = taun * n * T / 2
    return (dlnC_dlnT = dlnC_dlnT,
            ηbar = ηbar,                          # tau_n P0/2, P0 = nT
            ηbarα = ηbar,                         # the D alpha coefficient, eta_bar d_alpha ln I31
            λa   = Ds * (mq^2 + 6 * T * h) / (2 * T),
            Aco  = mq^2 + 5 * T * h,              # A = I41/P0
            ABr  = 5 + z * K2 / K3)               # A/B = dln I31/dln T
end

"""
    hq2d_m2_coeffs(T, n, taun, Ds, h, mq, mix) -> NamedTuple

The same block for a MULTI-SPECIES charm list, `mix = hq_m2_mixture(T, α, list, mq)`
(EquationOfState/heavyquark.jl): A/B, d ln C/d ln T and the m² of A from the list's
species sums, and the D α coefficient η̄ ∂_α ln I31, as in the two-way system. A
one-species list of mass `mq` uses the single-quark method above.
"""
function hq2d_m2_coeffs(T, n, taun, Ds, h, mq, mix)
    ηbar = taun * n * T / 2
    return (dlnC_dlnT = mix.dlnC_dlnT,
            ηbar = ηbar,
            ηbarα = ηbar * mix.aI31,
            λa   = Ds * (mix.m2bar + 6 * T * h) / (2 * T),
            Aco  = mix.m2bar + 5 * T * h,
            ABr  = mix.ABr)
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

so what must be subtracted from the LHS source is `tau_M (u^i c^j + u^j c^i)`.
Returned here as that quantity, already symmetrised; the caller subtracts it.

A 1-D reduction that maps the 2-D rate onto `p_l` with the projection held fixed
(`l_mu l_nu d_tau pi^{mu nu}`) misses this term: `l` depends on `u(tau)`, and
`d_tau(l l pi) - l l d_tau pi` is this term.

Derived symbolically: the orthogonality identity, the closed form on the
transverse block, and the absence of a trace remainder. Size: 1.6e-4 to 3.4e-3 of the rate.
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

⚠ `omega` is not `d_x u^y - d_y u^x`. A boost-invariant transverse flow carries
the acceleration terms too: `omega^{xy} = (1/2)(d_xu^y - d_yu^x + u^x a^y - u^y a^x)`.

Every eta-component of omega vanishes by boost invariance, so only the (tau,x,y)
block enters, and

    2 pi^{l<i} omega_l^{j>} = (pi g omega) + (pi g omega)^T

is ALREADY symmetric, traceless and u-orthogonal (pi symmetric, omega
antisymmetric, both u-orthogonal), so no projection is needed. All four
properties are checked symbolically.

⚠ Off by default (`vorticity = false`), so that existing results are unchanged.
It vanishes identically in 1+1D, since a radial flow has no transverse
vorticity, so no 1-D limit can test it. Its sign and normalisation are fixed by
the identity sigma-coupling + omega-coupling = the full contraction
`pi^{l<i} NablaPerp_l u^{j>}`, which ties it to the sigma coupling.

On a realistic Pb+Pb event `|omega|/|sigma|` has a median of ~2.8e-2 (converged
in resolution), and the term moves pi_Q by ~1e-2.
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
drags in terms carrying gradients of U. Derived symbolically; it reproduces the
1-D `(grad nu)_ll` to 1e-15.
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
                 dxa, dya, dta, dtnx, dxnx, dynx, dtny, dxny, dyny,
                 n, taun, Ds, h, hp, tauM, etaM, mq;
                 projector = true, vorticity = false,
                 cf = hq2d_m2_coeffs(T, n, taun, Ds, h, mq))
        -> (dpxx, dpxy, dpyy, dPiQ)

`d_tau` of the transverse second-moment block and the trace, under the CONSISTENT
closure. Returned as explicit rates rather than matrix rows because rows 4-8
relax diagonally (`tau_M u^tau d_tau X + src = 0`) and carry no derivative
couplings among themselves — the same structure the 1-D `hqc_m2_rows` has, where
`at3[i, 2+i] = tau_M u^tau` and the moment block is diagonal.

`piQ^eta_eta` is NOT returned: it is redundant. The driver keeps it on the
tracelessness value `hq2d_piQ_eta(...)` by a stiff relaxation in the RHS (its
row 7, rate 100/fm) — there is no post-step projection.

`cf` is the scalar coefficient block; the driver computes it once per cell and
passes it to all ten of its calls (it costs three Bessel functions).

`X = (alpha, nu^x, nu^y, pxx, pxy, pyy, peta, PiQ)`.

Two terms are switchable (see their helpers):

  * `projector = true` (default on): the Delta-projection of the comoving
    derivative, `hq2d_proj_corr`. It belongs in the equation; turning it off
    reproduces earlier results.
  * `vorticity = false` (default off): the DNMR corotation term,
    `hq2d_vorticity`. It also belongs in the equation; it is off so that
    existing results are unchanged.
"""
function hq2d_m2_rows(X, tau, T, ux, uy, dtT, dxT, dyT,
                      dtux, dxux, dyux, dtuy, dxuy, dyuy,
                      dxa, dya, dta, dtnx, dxnx, dynx, dtny, dxny, dyny,
                      n, taun, Ds, h, hp, tauM, etaM, mq;
                      projector::Bool = true, vorticity::Bool = false,
                      cf = hq2d_m2_coeffs(T, n, taun, Ds, h, mq))
    uT = sqrt(1 + ux^2 + uy^2)
    nux = X[2]; nuy = X[3]
    pxx = X[4]; pxy = X[5]; pyy = X[6]; pet = X[7]; PiQ = X[8]
    tauM > 0 || return (zero(T), zero(T), zero(T), zero(T))

    ηQ = etaM
    ζQ = (5//3) * ηQ

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
          (Ds*cf.Aco/(3*T))*adn + (5//6)*(Ds/T)*ngT +
          (hasproperty(cf, :ηbarα) ? cf.ηbarα : cf.ηbar)*Dα

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

⚠ Only c_M = 0 is supported, and rows 4-8 are then passive.
The second moments read the current and the background but feed back into
neither, so the 8×8 is block triangular: rows 1-3 are a 3-field first moment,
unchanged, and rows 4-8 hang off them. Rows 1-3 are `hq2d_matrices_consistent`
(the full-grad-P row, the 2-D analogue of `hqc_matrices`) by default, and
`hq2d_matrices` (the grad-alpha-only row of `matrix2d_HQ_BG!`) with
`consistent_fm = false`. It is the same structure the 1-D
`matrix1d_visc_HQ_BG_consistent_m2_5f!` has at its default c_M = 0.

τ_n, τ_M, η_M, h and D_s come from `hq_transport`, so `HQ_TAUN_SCALE` and
`FLUIDUM_HQ_MOMENTUM_DIMS` act on this driver exactly as on the 1-D ones.
"""
function matrix2d_HQ_BG_m2!(; dmn_eps = 1e-6, background_fields, α_max = 200.0,
                              consistent_fm::Bool = true,
                              projector::Bool = true, vorticity::Bool = false,
                              mixture::Bool = HQ_M2_MIXTURE[])
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_HQ_BG_m2!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields,
        α_max = α_max, consistent_fm = consistent_fm,
        projector = projector, vorticity = vorticity, mixture = mixture)
end

# `mixture` (default `HQ_M2_MIXTURE[]`, env FLUIDUM_HQ_M2_MIXTURE, true): for a multi-species charm list
# τ_M, d ln C/d ln T, A/B and A come from the list's species sums (`hq_m2_mixture`), as in the two-way
# system. false, or a one-species list of mass diffusion.mass: the single-quark forms.
function matrix2d_HQ_BG_m2!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0,
        consistent_fm::Bool = true,
        projector::Bool = true, vorticity::Bool = false,
        mixture::Bool = HQ_M2_MIXTURE[])

    T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy =
        background_fields(tau, X[1], X[2])

    α_safe = clamp(ϕ[1], -α_max, α_max)
    uT = sqrt(1 + ux^2 + uy^2)

    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_da = thermo.pressure_derivative
    dn_da += dmn_eps

    κ  = diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    # τ_n, τ_M, η_M, h, D_s: one τ_n for both moments, so τ_n = D_s h/T holds
    # under every diagnostic knob.
    tr   = hq_transport(T, α_safe, params.eos, params.diffusion)
    taun = tr.taun; tauM0 = tr.tauM; etaM = tr.etaM; h = tr.h; Ds = tr.Ds
    hp   = hqc_dh_dT(T, α_safe, params.eos, params.diffusion)

    # ── rows 1-3: the first moment ─────────────────────────────────────────────
    # `hq2d_matrices_consistent` (default) is the full-grad-P row,
    #     (D_s/T) Delta^mu_lambda Nabla_nu T_Q^{nu lambda} + nu^mu = 0,
    # generated by deriv/derive_hq_consistent_2p1d.wls (its cylindrical reduction
    # reproduces the 1-D consistent source, and its at-rest limit matches the closed
    # form to 2.2e-16). `consistent_fm = false` gives the grad-alpha-only row, whose
    # temperature-gradient channel is 0; on a realistic Pb+Pb cell that channel is
    # 1.8x the alpha drive.
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
    # The scalar coefficients are the same for all ten calls below: built once.
    # (tauM and cf assigned ONCE: the m2rate closure below captures them, and a captured variable
    # assigned on two paths is boxed — 17 kB per cell)
    tauM, cf = if mixture && !hq_m2_single(params.eos.hadron_list, params.diffusion.mass)
        mix = hq_m2_mixture(T, α_safe, params.eos.hadron_list, params.diffusion.mass)
        (tauM0 * mix.Cratio, hq2d_m2_coeffs(T, n, taun, Ds, h, params.diffusion.mass, mix))
    else
        (tauM0, hq2d_m2_coeffs(T, n, taun, Ds, h, params.diffusion.mass))
    end

    # ⚠ The alpha and nu gradients go in the matrix, not the source. Rows 4-8
    # depend on d_tau alpha, d_i alpha and d nu; in this convention those are
    # columns of At/Ax/Ay (columns 1, 2, 3), because alpha, nu^x and nu^y are
    # evolved fields of this system. Putting them in `src` with their current
    # values frozen would lag them by a step and hide a real coupling from the
    # characteristic analysis.
    #
    # `hq2d_m2_rows` returns explicit RATES, so the columns are read off by
    # finite-differencing it in each gradient. They are exactly linear in every
    # one of them (the system is quasi-linear), so a single difference is exact.
    m2rate(dta, dxa, dya, dtnx, dxnx, dynx, dtny, dxny, dyny) =
        hq2d_m2_rows((α_safe, ϕ[2], ϕ[3], ϕ[4], ϕ[5], ϕ[6], ϕ[7], ϕ[8]),
                     tau, T, ux, uy, dtT, dxT, dyT,
                     dtux, dxux, dyux, dtuy, dxuy, dyuy,
                     dxa, dya, dta, dtnx, dxnx, dynx, dtny, dxny, dyny,
                     n, taun, Ds, h, hp, tauM, etaM, params.diffusion.mass;
                     projector = projector, vorticity = vorticity, cf = cf)

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
    # relaxed toward its tracelessness value `hq2d_piQ_eta(...)`, which keeps it on
    # the constraint without adding a constraint equation. Tracelessness is
    # algebraic, not a slow mode: at a rate of 1/fm the field lags
    # visibly (piQ^eta/(pxx+pyy) = -0.34 instead of -1).
    #
    # `RELAX_RATE` sets how fast the constraint is enforced. It must be fast
    # compared with tau_M (the physical relaxation of the moments themselves), so
    # that row 7 tracks the transverse block instead of lagging it. 100/fm =
    # 1/(0.01 fm) is ~25x faster than a typical tau_M ~ 0.25 fm here and still
    # finite, which keeps the system an ODE rather than a DAE. It is a STIFF term
    # of the RHS, not a post-step projection: an explicit integrator's step is
    # stability-capped near ~0.03 fm by it (typical runs use dtmax = 0.01).
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
