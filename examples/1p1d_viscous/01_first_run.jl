#!/usr/bin/env julia
#=
01 — the minimal correct 1+1D viscous run, and the check to do before any other.

Everything here is the API you will copy: field layout, discretisation, fluid properties, initial
condition, driver, readout. The physics is deliberately trivial — a radially UNIFORM state at rest,
which must evolve as Bjorken — because that is the one configuration where 1+1D has a closed-form
answer, so it separates "the solver works" from "my initial condition is what I think it is".

This is the PRODUCTION kernel: `matrix1d_visc!` is what LangevinPaper1, AttractorHydro and
SpectraDiagnostic stand on, and the charm driver `matrix1d_visc_HQ!` (example 03) reuses its rows.

CONVENTIONS, which this file also pins:
  phi = (T, u^r, pi^phi_phi, pi^eta_eta, Pi)
  pi^eta_eta is the MIXED Milne component tau^2 pi^{eta eta} (NOT pi^{eta eta}), pi^phi_phi the
  mixed r^2 pi^{phi phi}.  Navier-Stokes fixed points on a Bjorken background:
      pi^eta_eta -> -4 eta/(3 tau),   pi^phi_phi -> +2 eta/(3 tau),   Pi -> -zeta/tau.
  The relaxation is ISRAEL-STEWART + delta_pipi pi theta, delta_pipi = (4/3) tau_pi (the DNMR value,
  default since 2026-10-13; FLUIDUM_DELTA_PIPI=0 restores plain IS); no tau_pipi, no lambda_piPi.
  The reference below follows the switch; the two delta_pipi values differ by ~20 % in pi^eta_eta on
  the state reached below, so the agreement pins which one the kernel carries.

THE AXIS IS NOT A BOUNDARY.  r = 0 is a symmetry axis: `ur` is declared ODD there and every other
field EVEN (`viscous_1d`).  Nothing in this file exercises that — a uniform state has u^r = 0
everywhere — but it is the first thing to get wrong on a real profile, so see example 02.
=#
ENV["GKSwstype"] = "100"
using Fluidum, Printf, OrdinaryDiffEq, Plots
# Figure style, shared by every example in BOTH packages (2026-09-14). dpi 200 and the larger
# fonts are for the GitHub READMEs: they render an image at container width (~900 px), so a
# 1700-px-wide figure is downscaled and 8 pt tick labels turn to mush. Keep the two packages
# identical — a reader comparing FiVo and Fluidum plots should not be reading two house styles.
gr(); default(; fontfamily = "sans-serif", framestyle = :box, grid = false, dpi = 200, lw = 2.2,
               titlefontsize = 11, guidefontsize = 10, tickfontsize = 9, legendfontsize = 8,
               foreground_color_legend = nothing, background_color_legend = RGBA(1,1,1,0.75),
               left_margin = 9Plots.mm, bottom_margin = 6Plots.mm, right_margin = 3Plots.mm,
               colorbar_titlefontsize = 8,   # ⚠ at 10 pt the colorbar TITLE overlaps its own
                                             # tick labels on every map panel (measured on ex10)
               top_margin = 2Plots.mm)
# ⚠ the margins are NOT cosmetic. Raising the font sizes above without them silently DROPS the axis
# labels and clips the y-label off the left edge — GR gives the axis whatever space is left after the
# panel, and at dpi 200 with 10 pt guides there is none. Measured on ex01: "τ [fm/c]" vanished from
# all three panels and reappeared only once the margins were set.
const F   = Fluidum
const FIG = joinpath(@__DIR__, "figures"); isdir(FIG) || mkpath(FIG)

const EOS  = FluiduMEoS()
const PAR  = F.FluidProperties(EOS, QGPViscosity(0.1, 0.2),        # eta/s = 0.1, tau_pi = 5(eta/s)/T
                                    SimpleBulkViscosity(0.083, 15.0),
                                    ZeroDiffusion())
const T0, TAU0, TAUF = 0.45, 0.4, 8.0

ent(T) = F.pressure_derivative(T, Val(1), EOS)
dst(T) = F.pressure_derivative(T, Val(2), EOS)
ene(T) = T*ent(T) - F.pressure(T, EOS)

# ── the run ─────────────────────────────────────────────────────────────────────────────────────
# `viscous_1d()` is the 5-field layout; `OriginInterval(N, R)` is a radial grid whose first cell
# centre sits at R/(2N), i.e. the axis is BETWEEN cells and is never evaluated directly.
disc = CartesianDiscretization(F.OriginInterval(24, 20.0))
d    = DiscreteFields(F.viscous_1d(), disc, Float64)
phi  = set_array(r -> T0, :temperature, d)
for f in (:ur, :piphiphi, :pietaeta, :piB); set_array!(phi, r -> 0.0, f, d); end
# ⚠ `set_array` returns a NEW array for the first field and `set_array!` writes into it. Building
# the IC with `set_array` twice silently discards the first field.
sol = F.oneshoot(d, F.matrix1d_visc!, PAR, phi, (TAU0, TAUF); reltol = 1e-10, abstol = 1e-12)

# ── the reference: the 0+1D Israel-Stewart system, written from the physics ──────────────────────
function rhs!(du, U, p, tau)
    T, pieta, Pi = U
    dpt  = ent(T); dptt = dst(T)
    eta  = F.viscosity(T, dpt, PAR.shear);       taupi = F.τ_shear(T, dpt, PAR.shear)
    zet  = F.bulk_viscosity(T, dpt, PAR.bulk);   tauPi = F.τ_bulk(T, dpt, dptt, PAR.bulk)
    du[1] = (-(ene(T) + F.pressure(T, EOS) + Pi + pieta)/tau)/(T*dst(T))   # de/dT = T d2P/dT2
    # theta = 1/tau on Bjorken; delta_pipi = Fluidum.DELTA_PIPI[] * (4/3) tau_pi, as the kernel runs
    du[2] = (-pieta - 4eta/(3tau) - F.DELTA_PIPI[]*(4/3)*taupi*pieta/tau)/taupi
    du[3] = (-Pi    - zet/tau)/tauPi
end
ref = solve(ODEProblem(rhs!, [T0, 0.0, 0.0], (TAU0, TAUF)), Tsit5(); reltol = 1e-12, abstol = 1e-14)

# ── the readout ─────────────────────────────────────────────────────────────────────────────────
# `sol` is an ODESolution over the whole field array: sol(tau)[ifield, icell]. Take a mid-grid cell
# so neither the axis nor the outer ghost can contribute.
k   = size(sol[end], 2) ÷ 2
taus = range(TAU0, TAUF; length = 200)
Tn, pn, Bn = ([sol(t)[i, k] for t in taus] for i in (1, 4, 5))
Tr, pr, Br = ([ref(t)[i]    for t in taus] for i in (1, 2, 3))

@printf("\n  at tau = %.1f fm/c, mid-grid cell\n", TAUF)
@printf("    T          solver %.8f   reference %.8f   rel %.2e\n", Tn[end], Tr[end], abs(Tn[end]-Tr[end])/Tr[end])
@printf("    pi^eta_eta solver %.8e   reference %.8e   rel %.2e\n", pn[end], pr[end], abs(pn[end]-pr[end])/abs(pr[end]))
@printf("    Pi         solver %.8e   reference %.8e   rel %.2e\n", Bn[end], Br[end], abs(Bn[end]-Br[end])/abs(Br[end]))
@printf("    |pi^eta_eta| / e = %.3f at this time; it peaks near tau0, where the gradients are\n",
        abs(pn[end])/ene(Tn[end]))
@printf("    |pi^eta_eta| / e = %.3f at tau = %.1f — that is the amplitude the agreement above is measured at\n",
        abs(sol(TAU0 + 0.2)[4, k])/ene(sol(TAU0 + 0.2)[1, k]), TAU0 + 0.2)
# pi^phi_phi = -pi^eta_eta/2 is NOT imposed: it comes out of two independently evolved rows.
@printf("    pi^phi_phi + pi^eta_eta/2 = %.2e  (relative %.2e)\n",
        sol[end][3, k] + pn[end]/2, abs(sol[end][3, k] + pn[end]/2)/abs(pn[end]))
# uniformity and rest are exact fixed points of the discretisation
@printf("    T spread across the grid %.2e,  max|u^r| %.2e\n",
        maximum(abs, sol[end][1, :] .- Tn[end])/Tn[end], maximum(abs, sol[end][2, :]))

plt = plot(layout = (1, 3), size = (1150, 440))
plot!(plt[1], taus, Tn; label = "solver", xlabel = "τ [fm/c]", ylabel = "T [GeV]", title = "temperature")
plot!(plt[1], taus, Tr; label = "0+1D MIS", ls = :dash)
plot!(plt[2], taus, pn; label = "solver", xlabel = "τ [fm/c]", ylabel = "π^η_η [GeV/fm³]", title = "shear")
plot!(plt[2], taus, pr; label = "0+1D MIS", ls = :dash)
plot!(plt[3], taus, Bn; label = "solver", xlabel = "τ [fm/c]", ylabel = "Π [GeV/fm³]", title = "bulk")
plot!(plt[3], taus, Br; label = "0+1D MIS", ls = :dash)
savefig(plt, joinpath(FIG, "ex01_bjorken.png"))
println("\n  wrote ", joinpath(FIG, "ex01_bjorken.png"))
