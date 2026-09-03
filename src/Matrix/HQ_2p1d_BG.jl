# ==============================================================================
# HQ_2p1d_BG.jl — THE CHARM SECTOR ON A PRESCRIBED 2+1D BACKGROUND.
#
# The 2-D twin of HQ_const_BG.jl.  Fields phi = (alpha, nu^x, nu^y) on Milne
# (tau, x, y, eta_s), boost invariant, evolved on a background T(tau,x,y) and
# u^i(tau,x,y) handed in as a closure — exactly the structure every production
# driver on this branch uses in 1-D (LP1, AttractorHydro, SpectraDiagnostic).
#
# The matrices come from Julia/tools/derive_hq_2p1d.wls and live in
# HQ_2p1d_BG_generated.jl.  That generator earns its trust by reproducing
# Fluidum's shipped 1-D cylindrical `one_d_viscous_matrix_fugacity_BG_only`
# ENTRY BY ENTRY from the same routine (its gate G4) before it is asked for the
# Cartesian chart.
#
# BACKGROUND CONTRACT.  `background_fields(tau, x, y)` returns the 12-tuple
#
#     (T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy)
#
# — the 2-D analogue of the 1-D 6-tuple (T, ur, dtT, drT, drur, dtur).  Note the
# 1-D tuple orders the u-derivatives (drur, dtur); this one groups by component,
# which is the only sane way to keep twelve of them straight.
#
# SCOPE.  Passive tracer: the charm does not feed T^{mu nu}.  Shipped
# (grad-alpha-only) first moment, i.e. the same closure as :HQ_const_BG — NOT the
# thermodynamically consistent full-grad-P one of HQ_const_BG_consistent.jl.
# ==============================================================================

"""
    matrix2d_HQ_BG!(; dmn_eps = 1e-6, background_fields, α_max = 200.0)

Closure constructor, matching `matrix1d_visc_HQ_BG!`.  Returns the
`(A_i, Source, ϕ, tau, X, params)` callable `basicupwinding` expects, with
`X = (x, y)` and `A_i` a 2-tuple of 3×3 matrices.
"""
function matrix2d_HQ_BG!(; dmn_eps = 1e-6, background_fields, α_max = 200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_HQ_BG!(
        A_i, Source, ϕ, tau, X, params;
        dmn_eps = dmn_eps, background_fields = background_fields, α_max = α_max)
end

function matrix2d_HQ_BG!(A_i, Source, ϕ, tau, X, params;
        dmn_eps = 1e-6, background_fields = nothing, α_max = 200.0)

    T, ux, uy, dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy =
        background_fields(tau, X[1], X[2])

    α_safe = clamp(ϕ[1], -α_max, α_max)

    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure                       # the charm-quark density (see charm_HRG.jl)
    dn_dT, dn_da = thermo.pressure_derivative
    dn_da += dmn_eps                          # the same regulator the 1-D matrices carry

    κ    = diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    taun = τ_diffusion_hadron(T, α_safe, params.eos, params.diffusion) * HQ_TAUN_SCALE[]

    At, Ax, Ay, src = hq2d_matrices(α_safe, ϕ[2], ϕ[3], tau, T, ux, uy,
                                    dtT, dxT, dyT, dtux, dxux, dyux, dtuy, dxuy, dyuy,
                                    n, dn_dT, dn_da, taun, κ)

    Ainv = inv(At)
    if !all(isfinite, Ainv)
        @warn "matrix2d_HQ_BG!: non-finite inv(A_t)" tau x=X[1] y=X[2] T ux uy n dn_da taun κ det_At=det(At) maxlog=5
    end

    A_mul_B!(A_i[1], Ainv, Ax)
    A_mul_B!(A_i[2], Ainv, Ay)
    jgemvavx!(Source, Ainv, src)
end

"""
    axisymmetric_bg_2d(bg1d)

Lift a 1-D cylindrical background closure `bg1d(tau, r) -> (T, ur, dtT, drT, drur, dtur)`
to the 2-D contract by rotating it into (x, y).  With `u^r` radial and everything
azimuthally symmetric,

    u^x = u^r x/r,  u^y = u^r y/r,
    d_x T = T'(r) x/r,   d_x u^x = u^r'(r) x²/r² + u^r y²/r³,   etc.

This is what makes the 2-D solver testable: run it on `axisymmetric_bg_2d(bg)` with an
axisymmetric initial condition and it must reproduce the 1-D cylindrical solve, which is
gated independently (Projects/FluidumValidation, gates B1–B3, S12).
"""
function axisymmetric_bg_2d(bg1d; rmin = 1e-9)
    function bg(tau, x, y)
        r  = max(sqrt(x*x + y*y), rmin)
        cx, cy = x/r, y/r
        T, ur, dtT, drT, drur, dtur = bg1d(tau, r)
        # scalars: gradient is radial
        dxT, dyT = drT*cx, drT*cy
        # vector u^i = u^r(r) x^i/r  ⇒  d_j u^i = (u^r' - u^r/r) x^i x_j /r² + (u^r/r) δ^i_j
        s = drur - ur/r
        dxux = s*cx*cx + ur/r
        dyux = s*cx*cy
        dxuy = s*cy*cx
        dyuy = s*cy*cy + ur/r
        (T, ur*cx, ur*cy, dtT, dxT, dyT, dtur*cx, dxux, dyux, dtur*cy, dxuy, dyuy)
    end
    return bg
end

# ==============================================================================
# THE 10-FIELD LIVE-BACKGROUND 2+1D SOLVER.
#
#   phi = (T, u^x, u^y, pi^yy, pi^zz, pi^xy, Pi_B, alpha, nu^x, nu^y)
#
# Rows 1-7 are the EXISTING 2-D viscous hydro (`two_d_viscous_matrix`, the block
# `matrix2d_visc!` already evolves).  Rows 8-10 are the charm equations of
# `hq2d_live_rows`, i.e. the same three equations as the background mode with the
# medium promoted from a closure to fields (gates G6a/G6b of derive_hq_2p1d.wls).
#
# PASSIVE TRACER (the house convention, and the operator's choice for this build):
# the charm does not feed T^{mu nu}, so the block (rows 1-7, cols 8-10) is ZERO.
# The reverse block (rows 8-10, cols 1-3) is NOT zero — the charm feels the medium.
#
# This REPLACES the old `matrix2d_visc_HQ!`, which cannot be used: its wrapper reads
# `phi[6]` as the fugacity while its own matrix documents `u[8]` as alpha, and with
# the documented layout its characteristic speeds reach 98-130 c with complex pairs
# (measured, Projects/FluidumValidation). The old entry point is left in place,
# unreferenced, rather than deleted — this repo keeps wrong turns, dated.
#
# Transport coefficients are taken through the SCALAR-entropy calls
# (`viscosity(T, dpt, shear)`) the 1-D `matrix1d_visc_HQ!` uses, not the
# `Thermodynamic`-object calls of `matrix2d_visc!` — the latter MethodError with
# `QGPViscosity`, the shear model every production run uses.
# ==============================================================================

"""
    matrix2d_visc_HQ_BG!(; dmn_eps = 1e-6, α_max = 200.0)

The 10-field 2+1D viscous hydro + passive charm system.  `X = (x, y)`; `A_i` is a
2-tuple of 10×10 matrices.
"""
function matrix2d_visc_HQ_BG!(; dmn_eps = 1e-6, α_max = 200.0)
    (A_i, Source, ϕ, tau, X, params) -> matrix2d_visc_HQ_BG!(
        A_i, Source, ϕ, tau, X, params; dmn_eps = dmn_eps, α_max = α_max)
end

const _HQ2D_COL = (1, 2, 3, 8, 9, 10)   # where hq2d_live_rows' six columns land

function matrix2d_visc_HQ_BG!(A_i, Source, ϕ, tau, X, params; dmn_eps = 1e-6, α_max = 200.0)
    T  = ϕ[1]
    dpt  = pressure_derivative(T, Val(1), params.eos)     # entropy
    dptt = pressure_derivative(T, Val(2), params.eos)
    η  = viscosity(T, dpt, params.shear)
    τS = τ_shear(T, dpt, params.shear)
    ζ  = bulk_viscosity(T, dpt, params.bulk)
    τB = τ_bulk(T, dpt, dptt, params.bulk)

    # rows 1-7: the medium, exactly as matrix2d_visc! evolves it
    uh = SVector{7}(ϕ[1], ϕ[2], ϕ[3], ϕ[4], ϕ[5], ϕ[6], ϕ[7])
    # `_active`, not `two_d_viscous_matrix`: this block inherits the ten wrong entries of the
    # shipped 2-D matrix and the `FLUIDUM_2D_DERIVED` switch that repairs them (see the note at
    # the top of Matrix/2d_viscous.jl).
    Ath, Axh, Ayh, srch = two_d_viscous_matrix_active(uh, tau, pressure(T, params.eos), dpt, dptt, ζ, η, τS, τB)

    # rows 8-10: the charm, reading the live medium
    α_safe = clamp(ϕ[8], -α_max, α_max)
    thermo = thermodynamic(T, α_safe, params.eos.hadron_list)
    n = thermo.pressure
    dn_dT, dn_da = thermo.pressure_derivative
    dn_da += dmn_eps
    κ    = diffusion_hadron(T, α_safe, params.eos, params.diffusion)
    taun = τ_diffusion_hadron(T, α_safe, params.eos, params.diffusion) * HQ_TAUN_SCALE[]
    AtL, AxL, AyL, srcL = hq2d_live_rows(α_safe, ϕ[9], ϕ[10], tau, T, ϕ[2], ϕ[3],
                                         n, dn_dT, dn_da, taun, κ)

    At = MMatrix{10,10,Float64}(undef); Ax = MMatrix{10,10,Float64}(undef)
    Ay = MMatrix{10,10,Float64}(undef); src = MVector{10,Float64}(undef)
    @inbounds for j in 1:10, i in 1:10
        if i <= 7 && j <= 7
            At[i,j] = Ath[i,j]; Ax[i,j] = Axh[i,j]; Ay[i,j] = Ayh[i,j]
        else
            At[i,j] = 0.0; Ax[i,j] = 0.0; Ay[i,j] = 0.0      # passive: (1:7, 8:10) stays zero
        end
    end
    @inbounds for i in 1:7; src[i] = srch[i]; end
    @inbounds for k in 1:6
        j = _HQ2D_COL[k]
        for i in 1:3
            At[7+i, j] = AtL[i, k]; Ax[7+i, j] = AxL[i, k]; Ay[7+i, j] = AyL[i, k]
        end
    end
    @inbounds for i in 1:3; src[7+i] = srcL[i]; end

    Ainv = inv(SMatrix{10,10}(At))
    A_mul_B!(A_i[1], Ainv, SMatrix{10,10}(Ax))
    A_mul_B!(A_i[2], Ainv, SMatrix{10,10}(Ay))
    jgemvavx!(Source, Ainv, SVector{10}(src))
end
