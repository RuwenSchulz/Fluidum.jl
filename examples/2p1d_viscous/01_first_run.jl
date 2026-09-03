#!/usr/bin/env julia
#=
01 — the minimal correct 2+1D viscous run, and the check to do before any other.

Everything here is the API you will copy: field layout, discretisation, fluid properties, initial
condition, driver, readout. The physics is deliberately trivial — a transversely UNIFORM state,
which must evolve as Bjorken — because that is the one configuration where 2+1D has a closed-form
answer, so it separates "the solver works" from "my initial condition is what I think it is".

CONVENTIONS, which this file also pins:
  phi = (T, u^x, u^y, pi^{yy}, pi^{zz}, pi^{xy}, Pi)
  pi^{zz} is the MIXED Milne component pi^eta_eta = tau^2 pi^{eta eta}  (NOT pi^{eta eta}),
  pi^{yy} is the Cartesian pi^y_y.  Navier-Stokes fixed points on a Bjorken background:
      pi^eta_eta -> -4 eta/(3 tau),  pi^y_y -> +2 eta/(3 tau),  Pi -> -zeta/tau.
  The relaxation is PLAIN ISRAEL-STEWART: no delta_pipi, no tau_pipi, no lambda_piPi.  A
  delta_pipi = (4/3) tau_pi term of the kind FiVoHydro carries would move pi^eta_eta by ~20 % on
  the state reached below, and it is measurably absent.

USE THE CORRECTED MATRIX FOR NEW WORK: `FLUIDUM_2D_DERIVED=1`.  The shipped generated matrix has ten
wrong entries in its two shear rows (the pi^{xy} row is driven by the one-sided velocity gradient
instead of its symmetric part).  It is exact for irrotational transverse flow, which is why nothing
in THIS file changes -- Bjorken has no transverse flow at all -- but it is wrong wherever the flow
carries vorticity, and it makes the characteristic speeds anisotropic.  See the note at the top of
src/Matrix/2d_viscous.jl.
=#
ENV["GKSwstype"] = "100"
using Fluidum, Printf, OrdinaryDiffEq, Plots
gr(); default(; fontfamily = "sans-serif", framestyle = :box, grid = false, dpi = 150, lw = 2)
const F   = Fluidum
const FIG = joinpath(@__DIR__, "figures"); isdir(FIG) || mkpath(FIG)

const EOS  = FluiduMEoS()
const PAR  = F.FluidProperties(EOS, QGPViscosity(0.1, 0.2),          # eta/s = 0.1, tau_pi = 5(eta/s)/T
                                    SimpleBulkViscosity(0.083, 15.0),
                                    ZeroDiffusion())
const T0, TAU0, TAUF = 0.45, 0.4, 8.0

ent(T) = F.pressure_derivative(T, Val(1), EOS)
dst(T) = F.pressure_derivative(T, Val(2), EOS)
ene(T) = T*ent(T) - F.pressure(T, EOS)

# ── the run ─────────────────────────────────────────────────────────────────────────────────────
# `viscous_2d()` is the 7-field layout; every boundary is a plain ghost because a Cartesian
# transverse grid has no symmetry axis to impose a parity on.
disc = CartesianDiscretization(F.SymmetricInterval(32, 12.0), F.SymmetricInterval(32, 12.0))
d    = DiscreteFields(F.viscous_2d(), disc, Float64)
phi  = set_array((x, y) -> T0, :temperature, d)
for f in (:ux, :uy, :piyy, :pizz, :pixy, :piB); set_array!(phi, (x, y) -> 0.0, f, d); end
# ⚠ `set_array` returns a NEW array for the first field and `set_array!` writes into it. Building
# the IC with `set_array` twice silently discards the first field.
sol = F.oneshoot(d, F.matrix2d_visc!, PAR, phi, (TAU0, TAUF); reltol = 1e-10, abstol = 1e-12)

# ── the reference: the 0+1D Israel-Stewart system, written from the physics ──────────────────────
function rhs!(du, U, p, tau)
    T, pieta, Pi = U
    th  = thermodynamic(T, EOS)
    eta = F.viscosity(T, th, PAR.shear);      taupi = F.τ_shear(T, th, PAR.shear)
    zet = F.bulk_viscosity(T, th, PAR.bulk);  tauPi = F.τ_bulk(T, th, PAR.bulk)
    du[1] = (-(ene(T) + F.pressure(T, EOS) + Pi + pieta)/tau)/(T*dst(T))   # de/dT = T d2P/dT2
    du[2] = (-pieta - 4eta/(3tau))/taupi
    du[3] = (-Pi - zet/tau)/tauPi
end
ref = solve(ODEProblem(rhs!, [T0, 0.0, 0.0], (TAU0, TAUF)), Tsit5(); reltol = 1e-12, abstol = 1e-14)

u  = sol(TAUF); ic = 17
@printf("\n  at tau = %.1f fm/c      solver            reference         relative\n", TAUF)
for (k, lbl, m) in ((1, "T", 1), (5, "pi^eta_eta", 2), (7, "Pi", 3))
    @printf("    %-12s %16.10f  %16.10f  %9.2e\n", lbl, u[k,ic,ic], ref(TAUF)[m],
            abs(u[k,ic,ic] - ref(TAUF)[m])/abs(ref(TAUF)[m]))
end
@printf("    %-12s %16.10f  %16.10f  %9.2e   (tracelessness, DYNAMICAL: two separate rows)\n",
        "pi^y_y", u[4,ic,ic], -ref(TAUF)[2]/2, abs(u[4,ic,ic] + ref(TAUF)[2]/2)/abs(ref(TAUF)[2]/2))
sp = maximum(abs, u[1,2:33,2:33] .- u[1,ic,ic])/u[1,ic,ic]
mv = maximum(abs, u[2,2:33,2:33]) + maximum(abs, u[3,2:33,2:33])
@printf("\n  transverse spread of T  %.2e  (a uniform state must stay uniform)\n", sp)
@printf("  max |u|                 %.2e  (a uniform state must stay at rest)\n", mv)
@printf("  |delta_pipi| would have moved pi^eta_eta by ~%.0f %% here — it is absent\n",
        100*(4/3)*F.τ_shear(ref(TAUF)[1], thermodynamic(ref(TAUF)[1], EOS), PAR.shear)/TAUF)

taus = range(TAU0, TAUF; length = 200)
plt  = plot(layout = (1,3), size = (1150, 360))
for (j, (k, m, lbl)) in enumerate(((1,1,"T  [GeV]"), (5,2,"pi^eta_eta  [GeV/fm^3]"), (7,3,"Pi  [GeV/fm^3]")))
    plot!(plt[j], taus, [sol(t)[k,ic,ic] for t in taus]; label = "2+1D solver",
          xlabel = "tau  [fm/c]", ylabel = lbl)
    plot!(plt[j], taus, [ref(t)[m] for t in taus]; label = "0+1D IS reference", ls = :dash, lc = :black)
end
savefig(plt, joinpath(FIG, "ex01_bjorken.png"))
println("\n  -> ", joinpath(FIG, "ex01_bjorken.png"))
