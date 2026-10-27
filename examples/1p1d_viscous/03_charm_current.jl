#!/usr/bin/env julia
#=
03 — the charm sector: `matrix1d_visc_HQ!`, and the conservation law that checks it.

Examples 01 and 02 evolve the medium; this one adds the charm fugacity `alpha` and its diffusion
current `nu^r` on a Woods-Saxon fireball.

  phi = (T, u^r, pi^phi_phi, pi^eta_eta, Pi, alpha, nu^r)        `HQ_viscous_1d()`

Charm is a passive tracer. It does not feed T^{mu nu}: the (hydro, charm) block of the matrix is
zero, so the medium evolves as without charm and the transport problem is linear in the charm
density. The same background can be replayed under different charm closures.

`alpha = mu/T` is the fugacity; the charm density follows from it through the equation of state.
It is not bounded like a temperature: on a realistic initial condition it reaches +20 where the
density is large and runs to the floor in the dilute tail. Do not clamp T to a pedestal while
passing alpha unchanged: that multiplies n by n_eq(T_min)/n_eq(T_raw), which is e^95 at a
0.012 GeV floor.

Charm number. Nothing creates or destroys charm quarks on these time scales, so the integral of the
charm current over a Cauchy slice must be independent of the slice. `test_integral_cauchy` computes
it. The drift it reports is first-order discretisation: it halves under refinement. A leak would be
resolution-independent. For reference: -5.93 / -3.01 / -1.54 / -0.80 % at N_r = 150 / 300 / 600 /
1200, ratios 1.97, 1.95, 1.92.
=#
ENV["GKSwstype"] = "100"
using Fluidum, Printf, OrdinaryDiffEq, Plots
# Figure style: dpi 200 and larger fonts keep the figures legible when a README scales them down.
gr(); default(; fontfamily = "sans-serif", framestyle = :box, grid = false, dpi = 200, lw = 2.2,
               titlefontsize = 11, guidefontsize = 10, tickfontsize = 9, legendfontsize = 8,
               foreground_color_legend = nothing, background_color_legend = RGBA(1,1,1,0.75),
               left_margin = 9Plots.mm, bottom_margin = 6Plots.mm, right_margin = 3Plots.mm,
               colorbar_titlefontsize = 8,   # at 10 pt the colorbar title overlaps its tick labels
               top_margin = 2Plots.mm)
# The margins are needed: at these font sizes GR otherwise drops the axis labels and clips the
# y-label off the left edge.
const F   = Fluidum
const FIG = joinpath(@__DIR__, "figures"); isdir(FIG) || mkpath(FIG)

# The charm equation of state needs a resonance list and a ccbar normalisation; `HQdiffusion(D_s, m)`
# carries the diffusion coefficient and the charm mass.
const HQ  = Heavy_Quark(readresonancelist(), 30.0)
const PAR = F.FluidProperties(HQ, QGPViscosity(0.1, 0.2), SimpleBulkViscosity(0.083, 15.0),
                              HQdiffusion(0.1, 1.5))
const R, TAU0, TAUF = 20.0, 0.4, 8.0
# AWS = 2.0, not the 0.6 of example 02. The charm sector is far more sensitive to the sharpness of
# the fireball edge than the medium. With 0.6 the solve dies at tau ~ 3.3-4.2 with "dt was forced
# below floating point epsilon": alpha runs away (to about -729) where the charm density underflows,
# and the alpha row's time coefficient u^tau dn_da = u^tau n vanishes with it. `dmn_eps` keeps A_t
# invertible, and where n < dmn_eps the regulator is the coefficient. The 2-D solver fails the same
# way. The cause is the alpha runaway, not the regulator; a soft edge avoids it, it does not fix it.
const T0, TMIN, AWS, ALPHA0 = 0.45, 0.05, 2.0, -3.0
profile(r) = T0/(exp((r - 6.0)/AWS) + 1) + TMIN

function run(N; tauf = TAUF)
    disc = CartesianDiscretization(F.OriginInterval(N, R))
    d    = DiscreteFields(F.HQ_viscous_1d(), disc, Float64)
    phi  = set_array(profile, :temperature, d)
    for f in (:ur, :piphiphi, :pietaeta, :piB, :nur); set_array!(phi, r -> 0.0, f, d); end
    set_array!(phi, r -> ALPHA0, :α, d)
    sol  = F.oneshoot(d, F.matrix1d_visc_HQ!, PAR, phi, (TAU0, tauf); reltol = 1e-8, abstol = 1e-10)
    (disc, sol)
end
interior(a) = a[2:end-1]

# ── the run, and what came out ──────────────────────────────────────────────────────────────────
disc, sol = run(300)
fin = sol[end]
rs  = [R*(i - 0.5)/300 for i in 1:300]
T, ur, al, nur = (interior(fin[i, :]) for i in (1, 2, 6, 7))

println("\n  A DIFFUSING CHARM CURRENT ON A COOLING FIREBALL  (N_r = 300, tau ", TAU0, " → ", TAUF, ")")
println("  " * "-"^80)
@printf("  T      : %.4f at the axis → %.4f at r = 10 fm\n", T[1], T[findfirst(>(10.0), rs)])
@printf("  u^r    : 0 at the axis → %.3f at its peak (r = %.1f fm) — drawn in example 02\n",
        maximum(ur), rs[argmax(ur)])
@printf("  alpha  : %.3f … %.3f   (started flat at %.1f)\n", minimum(al), maximum(al), ALPHA0)
@printf("  nu^r   : max %.4e at r = %.2f fm — zero at the axis by parity, zero where alpha is flat\n",
        maximum(abs, nur), rs[argmax(abs.(nur))])
@printf("  |nu^r|/n at its peak cell: %.4f\n",
        maximum(abs, nur)/max(F.thermodynamic(T[argmax(abs.(nur))], al[argmax(abs.(nur))], HQ.hadron_list).pressure, 1e-30))
println("\n  nu^r is BUILT, not initialised: it started identically zero and is driven by the gradient")
println("  of alpha that the cooling and the radial flow create. A current that stayed zero would")
println("  mean the drive is not wired; one that appeared at the axis would mean the parity is wrong.")

# ── the conservation law ────────────────────────────────────────────────────────────────────────
println("\n  CHARM NUMBER ON CAUCHY SLICES  (must be slice-independent; the drift is discretisation)")
println("  " * "-"^80)
@printf("  %6s %14s %16s %14s\n", "N_r", "Q(tau=0.5)", "max drift [%]", "ratio")
function drift_ladder()
    prev = NaN                    # in a function on purpose: `prev` assigned in a top-level `for`
    for N in (150, 300, 600)      #   is a new local every iteration (soft scope), and the ratio
        dsc, s = run(N)           #   column would print "—" forever.
        q = [c[1] for c in F.test_integral_cauchy(s, dsc.grid, (0.5, 8.0), HQ; dt = 0.5)]
        drift = 100*maximum(abs.(q .- q[1]))/q[1]
        @printf("  %6d %14.6f %16.3f %14s\n", N, q[1], drift,
                isnan(prev) ? "—" : @sprintf("%.2f", prev/drift))
        prev = drift
    end
end
drift_ladder()
println("\n  The drift HALVES with each refinement — first order, which is what `basicupwinding` gives.")
println("  A genuine leak would not move with the grid at all. That contrast is the whole check.")

# ── the figure ──────────────────────────────────────────────────────────────────────────────────
plt = plot(layout = (1, 3), size = (1400, 440))
# T alone: u^r (dimensionless, up to ~1.8) on the same axis would flatten T. The radial flow is
# drawn in example 02 and printed above.
plot!(plt[1], rs, T; label = "", c = 1, xlabel = "r [fm]", ylabel = "T [GeV]",
      title = "medium at τ = $TAUF", xlims = (0, 14))
hline!(plt[1], [0.1565]; c = :gray, ls = :dot, label = "T_fo")
plot!(plt[2], rs, al;  label = "", xlabel = "r [fm]", ylabel = "α = μ/T",
      title = "charm fugacity", xlims = (0, 14))
plot!(plt[3], rs, nur; label = "", xlabel = "r [fm]", ylabel = "ν^r [fm⁻³]",
      title = "diffusion current", xlims = (0, 14))
hline!(plt[3], [0.0]; c = :gray, ls = :dot, label = "")
savefig(plt, joinpath(FIG, "ex03_charm_current.png"))
println("\n  wrote ", joinpath(FIG, "ex03_charm_current.png"))
