#!/usr/bin/env julia
#=
02 — a real radial fireball: the sector ladder, the axis, and what resolution buys you.

Example 01 had no radial structure, so it could not see what happens at the axis, or how much of
the radial flow is physics rather than the scheme's own diffusivity. Both are measured here on a
Woods-Saxon profile.

1. The sector ladder. Ideal, +shear, +shear+bulk, from one initial condition, so the difference is
   the sector. Shear resists the transverse gradient and slows the radial push; bulk resists the
   expansion itself. Both are reported as a change in the flow the fireball builds.

2. The axis. `viscous_1d` declares `ur` odd at r = 0 and every other field even. A parity error
   does not crash; it puts a kink in the first few cells that looks like a physical feature. The
   measure is |u^r| in the first cell relative to its peak: it must go to zero, faster as the grid
   refines. A parity error leaves it finite and resolution-independent.

3. Convergence order. `basicupwinding` is first-order accurate, so the error falls like 1/N and the
   scheme carries a numerical diffusivity of order dx. The self-convergence ladder below measures the
   observed order. It should sit near 1; near 0.5 the run is not resolving its own initial condition.

Block-average the fine grid onto the coarse one; do not sub-sample it. Cell centres sit at
(i-1/2)dr, so a 2N grid's cells 2i-1 and 2i tile the N grid's cell i, and sub-sampling compares
different points and reports the profile's own slope as an error.
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

const EOS  = FluiduMEoS()
const SH   = QGPViscosity(0.1, 0.2)
const BU   = SimpleBulkViscosity(0.083, 15.0)
const IDEAL = F.FluidProperties(EOS, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())
const SHEAR = F.FluidProperties(EOS, SH,              ZeroBulkViscosity(), ZeroDiffusion())
const BOTH  = F.FluidProperties(EOS, SH,              BU,                  ZeroDiffusion())

const R, TAU0, TAUF = 20.0, 0.4, 5.0
const T0, TMIN, AWS = 0.45, 0.05, 0.6
profile(r) = T0/(exp((r - 6.0)/AWS) + 1) + TMIN

"one solve; returns (grid r, final field array, solution)"
function run(N, par; tauf = TAUF)
    disc = CartesianDiscretization(F.OriginInterval(N, R))
    d    = DiscreteFields(F.viscous_1d(), disc, Float64)
    phi  = set_array(profile, :temperature, d)
    for f in (:ur, :piphiphi, :pietaeta, :piB); set_array!(phi, r -> 0.0, f, d); end
    sol  = F.oneshoot(d, F.matrix1d_visc!, par, phi, (TAU0, tauf); reltol = 1e-8, abstol = 1e-10)
    rs   = [R*(i - 0.5)/N for i in 1:N]
    (rs, sol[end], sol)
end
interior(a) = a[2:end-1]                      # one ghost cell per side

# ── 1. the sector ladder ────────────────────────────────────────────────────────────────────────
println("\n  SECTOR LADDER  (N_r = 300, tau ", TAU0, " → ", TAUF, ")")
println("  " * "-"^78)
# quote |pi|/e over cells above freeze-out only. In the dilute tail e -> 0 while pi relaxes on a
# finite clock, so the ratio there is a property of the vacuum, not of the fluid, and it runs past 1.
const TFO = 0.1565
@printf("  %-18s %12s %12s %16s\n", "sector", "max u^r", "T(r=0)", "max|pi^eta_eta|/e")
lad = Dict{String,Any}()
for (name, par) in (("ideal", IDEAL), ("+ shear", SHEAR), ("+ shear + bulk", BOTH))
    rs, fin, _ = run(300, par)
    T = interior(fin[1, :]); ur = interior(fin[2, :]); pe = interior(fin[4, :])
    e = [T[i]*F.pressure_derivative(T[i], Val(1), EOS) - F.pressure(T[i], EOS) for i in eachindex(T)]
    hot = findall(>(TFO), T)
    lad[name] = (rs, T, ur)
    @printf("  %-18s %12.4f %12.5f %16.4f\n", name, maximum(ur), T[1],
            isempty(hot) ? NaN : maximum(abs.(pe[hot])./e[hot]))
end
@printf("\n  shear costs %.1f %% of the ideal radial flow; bulk a further %.1f %%.\n",
        100*(1 - maximum(lad["+ shear"][3])/maximum(lad["ideal"][3])),
        100*(1 - maximum(lad["+ shear + bulk"][3])/maximum(lad["+ shear"][3])))
println("  Both are RESISTANCE terms, so both reduce it — a sector that increased the flow would be")
println("  a sign error, and that is the cheapest check this ladder buys.")

# ── 2. the axis ─────────────────────────────────────────────────────────────────────────────────
println("\n  THE AXIS  (u^r must vanish there, faster as the grid refines)")
println("  " * "-"^78)
@printf("  %6s %14s %16s %14s\n", "N_r", "dr [fm]", "|u^r| first cell", "/ max|u^r|")
Ns = (100, 200, 400, 800)
axis = Float64[]
for N in Ns
    _, fin, _ = run(N, BOTH)
    ur = interior(fin[2, :])
    push!(axis, abs(ur[1])/maximum(abs, ur))
    @printf("  %6d %14.4f %16.3e %14.3e\n", N, R/N, abs(ur[1]), axis[end])
end
@printf("\n  ratio falls %.1fx over an 8x refinement, i.e. like N^-%.2f — FASTER than the 1/N a bare\n",
        axis[1]/axis[end], log(axis[1]/axis[end])/log(8))
println("  u^r ~ r would give, because the peak moves outward with resolution as well. What matters is")
println("  that it falls MONOTONICALLY: an axis parity error leaves this flat in N, or growing.")

# ── 3. self-convergence ─────────────────────────────────────────────────────────────────────────
println("\n  SELF-CONVERGENCE OF T  (block-averaged onto the coarsest grid)")
println("  " * "-"^78)
coarsen(a, f) = [sum(@view a[(i-1)*f+1:i*f])/f for i in 1:length(a)÷f]
Tprof = Dict(N => interior(run(N, BOTH)[2][1, :]) for N in (100, 200, 400))
l2(a, b) = sqrt(sum((a .- b).^2)/length(a))
d12 = l2(Tprof[100], coarsen(Tprof[200], 2))
d23 = l2(Tprof[200], coarsen(Tprof[400], 2))
@printf("  ‖T(100) − T(200)‖ = %.4e\n  ‖T(200) − T(400)‖ = %.4e\n  observed order p = %.2f\n",
        d12, d23, log2(d12/d23))
println("  `basicupwinding` is first order, so p ≈ 1 is the right answer. p well below 1 means the")
println("  INITIAL CONDITION is still coming into resolution, not that the scheme is worse than it is.")

# ── the figure ──────────────────────────────────────────────────────────────────────────────────
plt = plot(layout = (1, 3), size = (1250, 440))
for (name, c) in (("ideal", 1), ("+ shear", 2), ("+ shear + bulk", 3))
    rs, T, ur = lad[name]
    plot!(plt[1], rs, T;  label = name, c = c, xlabel = "r [fm]", ylabel = "T [GeV]",
          title = "temperature at τ = $TAUF", xlims = (0, 12))
    plot!(plt[2], rs, ur; label = name, c = c, xlabel = "r [fm]", ylabel = "u^r",
          title = "radial flow", xlims = (0, 12))
end
plot!(plt[3], collect(Ns), axis; m = :o, xscale = :log10, yscale = :log10, label = "",
      xlabel = "N_r", ylabel = "|u^r|₍first cell₎ / max|u^r|", title = "the axis refines away")
savefig(plt, joinpath(FIG, "ex02_fireball_and_resolution.png"))
println("\n  wrote ", joinpath(FIG, "ex02_fireball_and_resolution.png"))
