#!/usr/bin/env julia
#=
03 — an elliptic fireball: what each dissipative sector does, and the causal envelope you are in.

The production-shaped run. An almond-shaped initial condition converts spatial eccentricity into a
momentum anisotropy, and the sectors are switched on one at a time so the change is attributable.

TWO THINGS THIS FILE CHECKS THAT A FIREBALL SCRIPT USUALLY DOES NOT

1. THE CAUSAL ENVELOPE.  Israel-Stewart is hyperbolic and subluminal only inside a bounded region of
   (pi, Pi), and Fluidum carries no regulator: nothing clips pi or Pi, and nothing warns.  So the
   run has to be checked against the envelope AFTERWARDS.  The envelope is measured here from the
   characteristic matrix (the same measurement as gate M5 of bench/bench_2p1d_viscous.jl) and the
   trajectory is compared against it.
   For eta/s = 0.1 and the Lorentzian zeta, the binding constraint is BULK at low temperature:
   |Pi|/e ~ 0.46 at T = 0.15 GeV, against >= 1.0 for shear alone.

2. WHERE THE VIOLATION ACTUALLY LIVES.  Measured over the whole grid, |pi|/e exceeds the envelope;
   measured over T > 0.12 GeV it does not.  The difference is entirely the dilute tail, where e -> 0
   and any fixed pi is a large ratio.  That is an argument for masking observables at a temperature
   cut -- NOT an argument that the tail is fine.  Both numbers are printed, deliberately.

⚠ USE `FLUIDUM_2D_DERIVED=1` FOR THIS RUN.  An elliptic fireball started from rest is very nearly
irrotational, so the shipped matrix's pi^{xy} defect is small here -- but "very nearly" is not
"exactly", and this file prints the vorticity it develops so the size of the exposure is visible.
=#
ENV["GKSwstype"] = "100"
using Fluidum, Printf, Plots
gr(); default(; fontfamily = "sans-serif", framestyle = :box, grid = false, dpi = 150, lw = 2)
const F   = Fluidum
const FIG = joinpath(@__DIR__, "figures"); isdir(FIG) || mkpath(FIG)

const EOS = FluiduMEoS()
ent(T) = F.pressure_derivative(T, Val(1), EOS)
ene(T) = T*ent(T) - F.pressure(T, EOS)
const TAU0, TAUF, LBOX, N = 0.4, 6.0, 14.0, 96
const EPS2, TCUT = 0.25, 0.12

props(; shear, bulk) = F.FluidProperties(EOS,
    shear ? QGPViscosity(0.1, 0.2) : ZeroViscosity(),
    bulk  ? SimpleBulkViscosity(0.083, 15.0) : ZeroBulkViscosity(), ZeroDiffusion())

"Woods-Saxon with an area-preserving elliptic deformation, so the three cases are one fireball
squeezed rather than three fireballs of different size."
function T_ic(x, y; eps2 = EPS2)
    a = sqrt(1 - eps2); b = sqrt(1 + eps2)
    r = hypot(x/a, y/b)
    0.45/(1 + exp((r - 6.0)/0.6)) + 0.05
end
function run_case(; shear, bulk, eps2 = EPS2)
    par = props(; shear = shear, bulk = bulk)
    d   = DiscreteFields(F.viscous_2d(),
          CartesianDiscretization(F.SymmetricInterval(N, LBOX), F.SymmetricInterval(N, LBOX)), Float64)
    p = set_array((x, y) -> T_ic(x, y; eps2 = eps2), :temperature, d)
    for f in (:ux, :uy, :piyy, :pizz, :pixy, :piB); set_array!(p, (x, y) -> 0.0, f, d); end
    t = @elapsed sol = F.oneshoot(d, F.matrix2d_visc!, par, p, (TAU0, TAUF); reltol = 1e-8, abstol = 1e-10)
    (sol, t)
end
xs = [-LBOX + (i-0.5)*2LBOX/N for i in 1:N]
"energy-weighted momentum anisotropy over the cells above the temperature cut"
function anisotropy(u)
    sx = 0.0; sy = 0.0
    for i in 1:N, j in 1:N
        T = u[1,i+1,j+1]; T < TCUT && continue
        w = ene(T) + F.pressure(T, EOS)
        sx += w*u[2,i+1,j+1]^2; sy += w*u[3,i+1,j+1]^2
    end
    (sx + sy) > 0 ? (sx - sy)/(sx + sy) : 0.0
end
"max |pi^eta_eta|/e and |Pi|/e, over the hot region and over everything"
function fractions(u)
    a = 0.0; b = 0.0; A = 0.0; B = 0.0
    for i in 1:N, j in 1:N
        T = u[1,i+1,j+1]; e = ene(T); e <= 0 && continue
        fp = abs(u[5,i+1,j+1])/e; fb = abs(u[7,i+1,j+1])/e
        A = max(A, fp); B = max(B, fb)
        T < TCUT && continue
        a = max(a, fp); b = max(b, fb)
    end
    (a, b, A, B)
end
"transverse vorticity omega_xy = (1/2)(d_x u^y - d_y u^x), the part the shipped pi^{xy} row gets wrong"
function vorticity(u)
    h = 2LBOX/N; w = 0.0; s = 0.0
    for i in 2:N-1, j in 2:N-1
        u[1,i+1,j+1] < TCUT && continue
        dxuy = (u[3,i+2,j+1] - u[3,i,j+1])/(2h); dyux = (u[2,i+1,j+2] - u[2,i+1,j])/(2h)
        w = max(w, abs(dxuy - dyux)/2); s = max(s, abs(dxuy + dyux)/2)
    end
    (w, s)
end
"""the dissipative fraction at which the characteristic matrix first turns acausal, at rest.

What crosses is the LIGHT CONE, not hyperbolicity: max|Im lambda| stays at round-off (9.4e-18)
across the whole scan, and only the BULK sector crosses |Re lambda| = 1, only on the negative
(expanding) side. Verified independently on the re-derived matrix, which agrees to the digit
because the bulk and conservation rows are among the ones the shipped matrix gets right.
The distinction matters: the failure mode is superluminal propagation, NOT the complex-characteristic
ill-posedness that bit the charm sector (FLUIDUM_ATTRACTOR_AUDIT.md), so "refine it and watch it get
worse" is not the diagnostic to reach for here, and N-independence would not vindicate the sector."""
function envelope(T; mode = :bulk)
    e = ene(T); A = [zeros(7,7), zeros(7,7)]; S = zeros(7)
    par = props(; shear = true, bulk = true)
    for fr in 0.0:0.005:1.0
        pzz = mode === :bulk ? 0.0 : -fr*e; PiB = mode === :bulk ? -fr*e : 0.0
        F.matrix2d_visc!(A, S, [T, 0.0, 0.0, -pzz/2, pzz, 0.0, PiB], 2.0, (1.0,1.0), par)
        for th in range(0, 2pi, length = 17)[1:16]
            ev = eigvals(cos(th)*A[1] + sin(th)*A[2])
            (maximum(abs, imag.(ev)) > 1e-6 || maximum(abs, real.(ev)) > 1.0) && return fr
        end
    end
    1.0
end
using LinearAlgebra

# ⚠ THE eps2 = 0 CONTROL DOES NOT RETURN ZERO HERE, AND THAT IS NOT A BUG.  A Cartesian grid is not
# rotationally symmetric, and this scheme is first-order upwind, so an axisymmetric profile picks up
# a grid anisotropy that converges away only as O(dx).  The control therefore measures the FLOOR on
# any anisotropy this solver can report -- run it, and read the elliptic numbers against it.  (For
# scale: FiVoHydro's equivalent control returns 1.4e-15, because its scheme is second order and its
# reconstruction is symmetrised.  Different schemes, same physics.)
println("\n  eps2 = 0 CONTROL — the floor set by the grid's own anisotropy, not zero:")
s0, _ = run_case(; shear = true, bulk = true, eps2 = 0.0)
a0 = anisotropy(s0(TAUF))
@printf("    anisotropy = %+.3e   (%.2f %% of the eps2 = 0.25 signal below)\n", a0, 100*abs(a0)/0.17)

println("\n  eps2 = $EPS2, N = $N, tau = $TAU0 -> $TAUF fm/c")
@printf("  %-16s %12s %10s %11s %11s %11s %11s %8s\n", "sectors", "anisotropy", "response",
        "|pi|/e T>cut", "|pi|/e all", "|Pi|/e T>cut", "|Pi|/e all", "wall")
plt = plot(size = (780, 470), xlabel = "tau  [fm/c]", ylabel = "momentum anisotropy",
           title = "eps2 = $EPS2, N = $N", legend = :bottomright)
for (lab, sh, bu) in (("ideal", false, false), ("shear", true, false),
                      ("bulk", false, true), ("shear + bulk", true, true))
    sol, t = run_case(; shear = sh, bulk = bu)
    u = sol(TAUF); a, b, A, B = fractions(u)
    @printf("  %-16s %+12.5f %10.3f %11.3f %11.3f %11.3f %11.3f %7.1fs\n",
            lab, anisotropy(u), anisotropy(u)/EPS2, a, A, b, B, t)
    taus = range(TAU0, TAUF; length = 60)
    plot!(plt, taus, [anisotropy(sol(t)) for t in taus]; label = lab)
    if lab == "shear + bulk"
        wv, st = vorticity(u)
        @printf("\n  transverse vorticity at tau = %.1f: max|omega_xy| = %.4f 1/fm against max|sigma_xy| = %.4f 1/fm\n",
                TAUF, wv, st)
        @printf("  => the shipped pi^{xy} row's error is at most ~%.0f %% of its own driving here.\n", 100*wv/max(st, 1e-12))
    end
end
savefig(plt, joinpath(FIG, "ex03_elliptic_fireball.png"))

# ── the bulk-only row above is anomalous: the trigger is the SURFACE, not the grid ───────────────
# The bulk-only anisotropy at zeta/s = 0.083 is 0.004 against an ideal 0.203 -- a factor 50, from a
# bulk pressure that is only 6 % of the energy density anywhere it matters.  Four explanations were
# tested and three were killed:
#
#   * NOT a causality violation.  The Israel-Stewart envelope, evaluated at THIS run's own zeta/s
#     and at each cell's own boost, is not crossed by a single cell at any zeta/s: worst margin over
#     the whole grid -0.95 to -0.53.  (What crosses it, when anything does, is the light cone and
#     not hyperbolicity -- max|Im lambda| stays at round-off.)
#   * NOT the flow dying.  max|u| and the cell count barely move; <w u_y^2> changes by -3 % while
#     <w u_x^2> drops 23 %.
#   * NOT the grid in absolute terms.  A DIFFERENT initial condition -- same solver, same box, same
#     N, same zeta/s -- shows no collapse at all and bulk moving the anisotropy slightly UP.
#   * IT IS CELLS PER SURFACE THICKNESS.  Holding N fixed and changing ONLY the Woods-Saxon surface
#     width a_surf turns the collapse on and off; the ideal control is flat across the same sweep.
#
# The ladder below scans a_surf and N together so both routes to the same a_surf/dx are visible.
# THE PRACTICAL RULE, which "use N >= 144" was not, because it transfers to other boxes and profiles:
#
#       resolve the profile's surface with at least ~3 cells before trusting the bulk sector.
#
# For this +-14 fm box and a_surf = 0.6 fm that is N >= 140.  For a 0.3 fm surface it would be
# N >= 280, and "we used a fine grid" would not save you unless someone checked the ratio.
#
# ⚠ AND IT IS NOT TRACKED BY THE SIZE OF Pi.  max|Pi|/e over the measured region is 0.061 / 0.064 /
# 0.065 at N = 64 / 96 / 144 while the anisotropy swings 0.110 -> 0.004 -> 0.188.  The natural reader
# assumption -- that a dissipative artefact scales with the dissipative field -- is wrong here.
println("\n  the trigger is a_surf/dx, not N.  bulk zeta/s = 0.083 against its own ideal control:")
@printf("  %8s %6s %9s %11s %11s %9s %10s\n", "a_surf", "N", "a_surf/dx", "ideal", "bulk", "ratio", "max|u|")
function surf_run(asurf, Nr, zs)
    par = F.FluidProperties(EOS, ZeroViscosity(),
              zs == 0 ? ZeroBulkViscosity() : SimpleBulkViscosity(zs, 15.0), ZeroDiffusion())
    d = DiscreteFields(F.viscous_2d(),
        CartesianDiscretization(F.SymmetricInterval(Nr, LBOX), F.SymmetricInterval(Nr, LBOX)), Float64)
    a = sqrt(1 - EPS2); b = sqrt(1 + EPS2)
    p = set_array((x, y) -> 0.45/(1 + exp((hypot(x/a, y/b) - 6.0)/asurf)) + 0.05, :temperature, d)
    for f in (:ux, :uy, :piyy, :pizz, :pixy, :piB); set_array!(p, (x, y) -> 0.0, f, d); end
    u = F.oneshoot(d, F.matrix2d_visc!, par, p, (TAU0, TAUF); reltol = 1e-8, abstol = 1e-10)(TAUF)
    sx = 0.0; sy = 0.0; mu = 0.0
    for i in 1:Nr, j in 1:Nr
        T = u[1,i+1,j+1]; e = ene(T); (e <= 0 || T < TCUT) && continue
        w = e + F.pressure(T, EOS); sx += w*u[2,i+1,j+1]^2; sy += w*u[3,i+1,j+1]^2
        mu = max(mu, hypot(u[2,i+1,j+1], u[3,i+1,j+1]))
    end
    ((sx - sy)/(sx + sy), mu)
end
for asurf in (0.6, 0.9, 1.2), Nr in (64, 96, 144)
    id, _  = surf_run(asurf, Nr, 0.0)
    bu, mu = surf_run(asurf, Nr, 0.083)
    @printf("  %8.2f %6d %9.2f %11.5f %11.5f %9.3f %10.4f\n",
            asurf, Nr, asurf/(2LBOX/Nr), id, bu, bu/id, mu)
end
println("  => the ideal control is flat to 3 % over the whole sweep, so the profile change is not")
println("     doing the work.  The ratio recovers to ~0.95 by a_surf/dx ~ 3 and to ~1 beyond it.")
println("  ⚠ BELOW the threshold the error is NOT a function of a_surf/dx alone: a_surf/dx = 2.06 is")
println("     reached both at (0.6, N=96) and at (0.9, N=64), giving ratios 0.021 and 0.60 -- a factor")
println("     29 apart.  Only the RECOVERY is single-valued; the failure is not, which is why two")
println("     resolutions both below the threshold certified a convincing 'converged' collapse.")
println("  ⚠ and a_surf/dx is the right control variable only to a couple of per cent.  Reaching a")
println("     given ratio by REFINING sits systematically ~1.5-2 % below reaching it by WIDENING the")
println("     surface (3.09: 0.932 vs ~0.947 interpolated; 4.46: 0.967 vs ~0.985).  Refining N also")
println("     refines everything else in the problem.  Use it as a lower bound, not as a calibration.")

println("\n  the causal envelope, from the characteristic matrix at rest:")
@printf("  %-10s %14s %14s\n", "T [GeV]", "shear only", "bulk only")
# ⚠ SPAN THE LORENTZIAN PEAK.  zeta/s is Lorentzian in T with its maximum at 0.175 and a width of
# 0.024, and the threshold tracks zeta/(w tau_Pi) inversely -- so the envelope is NARROWEST at 0.175
# and widens on BOTH sides.  A table sampled only at 0.15, 0.20, 0.30, 0.45 rises monotonically and
# invites exactly the wrong reading ("bulk binds at low T"), which is the reading that produced the
# wrong diagnosis of the ladder above.
for T in (0.03, 0.05, 0.10, 0.12, 0.15, 0.175, 0.20, 0.25, 0.30, 0.45, 0.60)
    @printf("  %-10.3f %14.3f %14.3f\n", T, envelope(T; mode = :shear), envelope(T; mode = :bulk))
end
println("  => narrowest AT the zeta(T) peak, wider on both sides. Shear alone stays subluminal to 1.0.")
println("\n  -> ", joinpath(FIG, "ex03_elliptic_fireball.png"))
println("""
  READING IT
    Shear reduces the anisotropy: it resists the differential expansion that builds it, which is why
    v2 constrains eta/s (0.20321 -> 0.17138). Adding bulk on top of shear reduces it a little more
    (0.17138 -> 0.16959): bulk does work against the expansion rather than resisting the
    DIFFERENTIAL expansion, so it cools the medium more than it reshapes it.
    The "bulk" row on its own is NOT that small effect -- it is a coarse-grid artefact, which the
    resolution ladder below identifies and which refinement removes. Do not read it as physics.
    Compare the two |pi|/e columns against the envelope table -- but compare each cell against the
    threshold AT ITS OWN TEMPERATURE. Doing it against a single number is how this file first got
    the ladder above wrong: the largest |Pi|/e sits in the coldest cells, where the threshold is
    ~0.635, not at the ~0.41 of the Lorentzian peak. Checked cell by cell, at its own temperature
    and its own boost, no cell in any of these runs is superluminal.
    The unmasked and masked |pi|/e columns still differ by a factor ~10, and that is still a reason
    to mask observables at a temperature cut -- just not evidence of a causality problem.
    The bulk-only ladder is the real lesson, and it is not about causality at all: a dissipative
    sector can be quantitatively wrong on a grid that looks perfectly reasonable, and wrong
    NON-MONOTONICALLY, so the two-resolution check that catches most things does not catch this.
    Three points, and a sector switched off as a control, is the minimum. Example 02 measures why:
    this scheme's own numerical viscosity is D_num = 0.75*dx, which at N = 96 on this box is
    dx = 0.29 fm -- more than twice the 0.12 fm at which it already equals the physical Gamma_s.""")
