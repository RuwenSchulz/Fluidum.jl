#!/usr/bin/env julia
# ════════════════════════════════════════════════════════════════════════════════════════════════
# bench/bench_1p1d_viscous.jl — THE 1+1D VISCOUS SOLVER AGAINST CLOSED FORMS.
#
# `matrix1d_visc!` (src/Matrix/1d_viscous_cilindrical.jl) evolves
#
#       φ = (T, u^r, π^φ_φ, π^η_η, Π)
#
# on a cylindrical, boost-invariant Milne slice. **This is the production kernel** — LangevinPaper1,
# AttractorHydro and SpectraDiagnostic all stand on the matrix this file checks, and its charm
# sibling `matrix1d_visc_HQ!` reuses the same five rows.
#
# WHY IT EXISTS (2026-09-14). The 2+1D solver got `bench/bench_2p1d_viscous.jl` on 2026-09-03 —
# 19 gates against closed forms — and the 1+1D one, which is the one every paper actually uses, had
# no equivalent inside the package at all. Its coverage lived entirely in
# `Julia/Projects/FluidumValidation/`, which `run_suite.jl` cannot discover, and in scattered
# asserts in `test/runtests.jl`. This file is the 1-D twin of the 2-D ladder: same structure, same
# conventions, every reference a closed form or an ODE integrated inside this file. Nothing reads
# the DataStore, so it is portable and deterministic.
#
#   julia --project=Julia Julia/Fluidum.jl/bench/bench_1p1d_viscous.jl
#   FLUIDUM_1DV_FULL=1 …    adds N = 800 to the axis-parity ladder in A
#
#   M  the characteristic structure, from the matrix alone — closed-form speeds, causal envelope
#   D  the DERIVATION control: `one_d_viscous_matrix_derived` reproduces the shipped kernel
#   B  BJORKEN   a radially uniform state must integrate the 0+1D Israel–Stewart equations EXACTLY
#   A  STRUCTURE the origin parity, and convergence of a real fireball
#   C  COST
#
# CONVENTIONS, pinned by B1 and used throughout
# ─────────────────────────────────────────────
#   `piphiphi` and `pietaeta` are MIXED components: π^η_η = τ²π^{ηη}, π^φ_φ = r²π^{φφ}. Same mixed
#   convention `viscous_2d` uses for `piyy`/`pizz`. Navier–Stokes fixed points on Bjorken:
#       π^η_η → −4η/(3τ),   π^φ_φ → +2η/(3τ),   Π → −ζ/τ.
#   The relaxation is PLAIN ISRAEL–STEWART — no δ_ππ, no τ_ππ, no λ_πΠ. B1 pins that to 1e-13; a
#   δ_ππ = (4/3)τ_π of the kind FiVoHydro carries would move π^η_η by ~20 % on the state reached.
# ════════════════════════════════════════════════════════════════════════════════════════════════
using Fluidum, LinearAlgebra, Printf, OrdinaryDiffEq
const F    = Fluidum
const FULL = get(ENV, "FLUIDUM_1DV_FULL", "0") == "1"

mutable struct Led; p::Int; f::Int; names::Vector{String}; end
const L = Led(0, 0, String[])
function gate(id, ok::Bool, msg)
    ok ? (L.p += 1) : (L.f += 1; push!(L.names, id))
    @printf("  [%s] %-5s %s\n", ok ? "PASS" : "FAIL", id, msg); ok
end
note(m)    = println("        ", m)
section(t) = (println(); println("── ", t, " ", "─"^max(0, 96 - length(t))))

const EOS   = FluiduMEoS()
const SHEAR = QGPViscosity(0.1, 0.2)                 # η/s = 0.1, τ_π = η/(s T C_s) = 5(η/s)/T̂
const BULK  = SimpleBulkViscosity(0.083, 15.0)
const PAR   = F.FluidProperties(EOS, SHEAR,           BULK,                ZeroDiffusion())
const PAR_I = F.FluidProperties(EOS, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())

prs(T)  = F.pressure(T, EOS)
ent(T)  = F.pressure_derivative(T, Val(1), EOS)
dst(T)  = F.pressure_derivative(T, Val(2), EOS)
ene(T)  = T*ent(T) - prs(T)
enth(T) = T*ent(T)
cs2(T)  = ent(T)/(T*dst(T))

"(η, τ_π, ζ, τ_Π) exactly as `matrix1d_visc!` builds them"
function coeffs(T, par = PAR)
    dpt = F.pressure_derivative(T, Val(1), par.eos)
    dptt = F.pressure_derivative(T, Val(2), par.eos)
    (F.viscosity(T, dpt, par.shear), F.τ_shear(T, dpt, par.shear),
     F.bulk_viscosity(T, dpt, par.bulk), F.τ_bulk(T, dpt, dptt, par.bulk))
end

"(A_t⁻¹A_r, A_t⁻¹s) at one state — what the solver actually integrates"
function sys(T, ur, pphi, peta, PiB, τ, r, par = PAR)
    A = [zeros(5,5)]; S = zeros(5)
    F.matrix1d_visc!(A, S, [T, ur, pphi, peta, PiB], τ, (r,), par)
    A[1], S
end

grid1d(N, R) = DiscreteFields(F.viscous_1d(),
                              CartesianDiscretization(F.OriginInterval(N, R)), Float64)

println("="^100)
println("Fluidum 1+1D VISCOUS HYDRO — closed-form validation ladder   (", FULL ? "FULL" : "fast", " mode)")
println("="^100)

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("M · the characteristic structure — closed forms, no PDE solve")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# In radial symmetry π^{μν} is diagonal in the parallel-transported triad (l, φ̂, η̂), so the two
# stored traceless channels fix the third: π^l_l = −(π^φ_φ + π^η_η). M1 checks that the tensor the
# solver's three numbers describe really is traceless and u-orthogonal — the 1-D statement of the
# 2-D bench's M1.
function pi_full(ur, pphi, peta)
    ut = sqrt(1 + ur^2)
    pll = -(pphi + peta)                     # tracelessness in the triad
    # l^μ = (u^r, u^τ, 0, 0) is the unit vector along r, orthogonal to u.  In the (τ, r, φ̂, η̂)
    # orthonormal frame π = diag over (l, φ̂, η̂) lifted back through l:
    ptt = pll*ur^2;  ptr = pll*ur*ut;  prr = pll*ut^2
    [ptt ptr 0.0 0.0; ptr prr 0.0 0.0; 0.0 0.0 pphi 0.0; 0.0 0.0 0.0 peta]
end
function gate_M1()
    g = Diagonal([-1.0, 1.0, 1.0, 1.0]); wtr = 0.0; wor = 0.0
    for ur in (-2.0, -0.5, 0.0, 0.7, 3.0), (pphi, peta) in ((0.1, -0.2), (-0.03, 0.07), (0.5, 0.5))
        Π = pi_full(ur, pphi, peta); u = [sqrt(1+ur^2), ur, 0.0, 0.0]
        sc = max(maximum(abs, Π), 1e-30)
        wtr = max(wtr, abs(tr(g*Π))/sc); wor = max(wor, maximum(abs, (u'*g*Π)')/sc)
    end
    gate("M1", wtr < 1e-13 && wor < 1e-13,
        @sprintf("the 2 stored channels close onto a traceless, u-orthogonal π^{μν} (|π^μ_μ| %.1e, |u_μπ^{μν}| %.1e, relative)", wtr, wor))
end
gate_M1()

# At rest and with π = Π = 0 the largest characteristic of Israel–Stewart is the MIS front
#     v_front = sqrt(c_s² + (4/3)η/(w τ_π) + ζ/(w τ_Π)),   w = e + p = T s.
function gate_M2()
    worst = 0.0; vs = Float64[]
    for T in (0.12, 0.2, 0.3, 0.45, 0.6)
        η, τπ, ζ, τΠ = coeffs(T); w = enth(T)
        Ar, _ = sys(T, 0.0, 0.0, 0.0, 0.0, 1.0, 5.0)
        λ = maximum(abs, real.(eigvals(Ar)))
        v = sqrt(cs2(T) + (4/3)*η/(w*τπ) + ζ/(w*τΠ))
        push!(vs, v); worst = max(worst, abs(λ - v)/v)
    end
    gate("M2", worst < 1e-12,
        @sprintf("at rest the largest characteristic ≡ √(c_s² + (4/3)η/(wτ_π) + ζ/(wτ_Π)), the MIS front (%.2e; v = %.3f…%.3f)",
                 worst, minimum(vs), maximum(vs)))
end
gate_M2()

# With every dissipative coefficient switched off the extreme characteristics must be the
# relativistic sound cone (v ± c_s)/(1 ± v c_s) — the ideal limit, boosted.
function gate_M3()
    worst = 0.0
    for T in (0.15, 0.3, 0.5), ur in (-1.5, -0.4, 0.0, 0.9, 2.5)
        Ar, _ = sys(T, ur, 0.0, 0.0, 0.0, 1.0, 5.0, PAR_I)
        λ = sort(real.(eigvals(Ar))); v = ur/sqrt(1 + ur^2); c = sqrt(cs2(T))
        lo = (v - c)/(1 - v*c); hi = (v + c)/(1 + v*c)
        worst = max(worst, abs(λ[1] - lo)/max(abs(lo), 1e-3), abs(λ[end] - hi)/max(abs(hi), 1e-3))
    end
    gate("M3", worst < 1e-11,
        @sprintf("with viscosity off the extreme characteristics ≡ the sound cone (v±c_s)/(1±v c_s) (%.2e)", worst))
end
gate_M3()

# A_t must stay invertible over the temperatures a fireball actually visits. Below that the EoS
# derivatives underflow and the system is singular — the 1-D half of §E1 of bench_fluidum_full.jl,
# which measured the floor at ≈0.02 GeV. Here we assert only the range production uses.
function gate_M4()
    worst = 0.0; Tw = 0.0
    for T in (0.05, 0.08, 0.12, 0.2, 0.3, 0.45, 0.6), ur in (0.0, 1.0)
        dpt = ent(T); dptt = dst(T); η, τπ, ζ, τΠ = coeffs(T)
        At, _, _ = F.one_d_viscous_matrix([T, ur, 0.0, 0.0, 0.0], 1.0, 5.0, prs(T), dpt, dptt, ζ, η, τπ, τΠ)
        c = cond(Matrix(At)); c > worst && (worst = c; Tw = T)
    end
    gate("M4", isfinite(worst) && worst < 1e8,
        @sprintf("A_t stays invertible over T = 0.05…0.6 GeV (worst cond %.2e at T = %.2f)", worst, Tw))
end
gate_M4()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("D · the derivation control — the shipped kernel against an independent re-derivation")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# `Julia/tools/derive_2p1d_viscous.wls` re-derives Israel–Stewart from scratch and emits
# `Matrix/viscous_generated.jl`.  Its 1-D routine is what makes the whole derivation trustworthy:
# it must reproduce the SHIPPED production kernel entry by entry.  That control has existed since
# 2026-09-03 but lived as gate G1 of a file named for the 2-D sector and was documented in no 1-D
# document, so nothing running the 1-D side ever exercised it. It runs here now.
function gate_D1()
    worst = 0.0; nst = 0
    for T in (0.12, 0.2, 0.35, 0.55), ur in (-1.3, -0.2, 0.0, 0.6, 2.2),
        (pphi, peta, PiB) in ((0.0,0.0,0.0), (0.02,-0.05,0.01), (-0.11,0.07,-0.03)),
        τ in (0.4, 2.0, 7.0), r in (0.05, 3.0, 12.0)
        η, τπ, ζ, τΠ = coeffs(T); φ = [T, ur, pphi, peta, PiB]
        a = F.one_d_viscous_matrix(        φ, τ, r, prs(T), ent(T), dst(T), ζ, η, τπ, τΠ)
        b = F.one_d_viscous_matrix_derived(φ, τ, r, prs(T), ent(T), dst(T), ζ, η, τπ, τΠ)
        for (x, y) in zip(a, b)
            sc = max(maximum(abs, x), 1e-30)
            worst = max(worst, maximum(abs, Array(x) .- Array(y))/sc)
        end
        nst += 1
    end
    gate("D1", worst < 1e-12,
        @sprintf("`one_d_viscous_matrix_derived` ≡ the shipped `one_d_viscous_matrix`, entry by entry over %d states (max rel %.2e)", nst, worst))
end
gate_D1()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("B · Bjorken — a radially uniform state IS the 0+1D Israel–Stewart system")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# The one configuration with a closed-form answer. A state with no radial dependence and u^r = 0
# must integrate
#   de/dτ = −(e + p + Π + π^η_η)/τ ;  τ_π Dπ^η_η = −π^η_η − 4η/(3τ) ;  τ_Π DΠ = −Π − ζ/τ
# and π^φ_φ = −π^η_η/2 must EMERGE from two independently evolved rows rather than be imposed.
function bjorken_reference(T0, τ0, τf)
    function rhs!(du, U, p, τ)
        T, peta, PiB = U
        η, τπ, ζ, τΠ = coeffs(T)
        du[1] = (-(ene(T) + prs(T) + PiB + peta)/τ)/(T*dst(T))       # de/dT = T d²p/dT²
        du[2] = (-peta - 4η/(3τ))/τπ
        du[3] = (-PiB  - ζ/τ)/τΠ
    end
    solve(ODEProblem(rhs!, [T0, 0.0, 0.0], (τ0, τf)), Tsit5(); reltol = 1e-12, abstol = 1e-14)
end

function gate_B()
    T0, τ0, τf = 0.45, 0.4, 6.0
    d   = grid1d(24, 20.0)
    phi = set_array(r -> T0, :temperature, d)
    for f in (:ur, :piphiphi, :pietaeta, :piB); set_array!(phi, r -> 0.0, f, d); end
    sol = F.oneshoot(d, F.matrix1d_visc!, PAR, phi, (τ0, τf); reltol = 1e-11, abstol = 1e-13)
    ref = bjorken_reference(T0, τ0, τf)

    fin = sol[end]; k = size(fin, 2) ÷ 2          # a mid-grid cell, away from both boundaries
    Tn, urn, pphin, petan, PiBn = (fin[i, k] for i in 1:5)
    Tr, petar, PiBr = ref(τf)
    eT  = abs(Tn - Tr)/Tr
    ep  = abs(petan - petar)/max(abs(petar), 1e-30)
    eB  = abs(PiBn - PiBr)/max(abs(PiBr), 1e-30)
    gate("B1", eT < 1e-10 && ep < 1e-9 && eB < 1e-9,
        @sprintf("a uniform state integrates the 0+1D MIS equations: T %.2e, π^η_η %.2e, Π %.2e (relative)", eT, ep, eB))
    note(@sprintf("state reached: T = %.5f, π^η_η = %.5e, Π = %.5e, |π^η_η|/e = %.3f",
                  Tn, petan, PiBn, abs(petan)/ene(Tn)))

    resid = abs(pphin + petan/2)/max(abs(petan), 1e-30)
    gate("B2", resid < 1e-10,
        @sprintf("π^φ_φ = −π^η_η/2 EMERGES from two independently evolved rows (residual %.2e)", resid))

    spread = maximum(abs, fin[1, :] .- fin[1, k])/Tn
    gate("B3", spread < 1e-12 && maximum(abs, fin[2, :]) < 1e-12,
        @sprintf("uniformity and rest are preserved to round-off (T spread %.1e, max|u^r| %.1e)",
                 spread, maximum(abs, fin[2, :])))

    # Ideal limit: with no dissipation τ·s is conserved on the Bjorken solution.
    phi2 = set_array(r -> T0, :temperature, d)
    for f in (:ur, :piphiphi, :pietaeta, :piB); set_array!(phi2, r -> 0.0, f, d); end
    s2  = F.oneshoot(d, F.matrix1d_visc!, PAR_I, phi2, (τ0, τf); reltol = 1e-11, abstol = 1e-13)
    e0  = τ0*ent(T0); e1 = τf*ent(s2[end][1, k])
    gate("B4", abs(e1 - e0)/e0 < 1e-10,
        @sprintf("with viscosity off, τ·s is conserved on the Bjorken solution (%.2e)", abs(e1-e0)/e0))
end
gate_B()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("A · structure — the origin, and a real fireball")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# r = 0 is a SYMMETRY AXIS, not a boundary: u^r is odd there and everything else even. A profile
# that is smooth and flat at the origin must stay so; a parity error shows as a kink in the first
# cells, which is exactly what this measures.
function fireball(N, R, τf; par = PAR, a = 0.6, T0 = 0.45, Tmin = 0.05)
    d   = grid1d(N, R)
    phi = set_array(r -> T0/(exp((r - 6.0)/a) + 1) + Tmin, :temperature, d)
    for f in (:ur, :piphiphi, :pietaeta, :piB); set_array!(phi, r -> 0.0, f, d); end
    (d, F.oneshoot(d, F.matrix1d_visc!, par, phi, (0.4, τf); reltol = 1e-8, abstol = 1e-10))
end

function gate_A()
    d, sol = fireball(200, 20.0, 4.0)
    fin = sol[end]
    gate("A1", all(isfinite, fin) && minimum(fin[1, :]) > 0,
        @sprintf("a Woods–Saxon fireball reaches τ = 4 with every field finite and T > 0 (T_max %.4f, max|u^r| %.4f)",
                 maximum(fin[1, :]), maximum(abs, fin[2, :])))

    # u^r must vanish at the axis and grow linearly out of it: |u^r| at the first cell should be
    # small compared with its peak, and shrink as the grid refines.
    firsts = Float64[]
    Ns = FULL ? (100, 200, 400, 800) : (100, 200, 400)
    for N in Ns
        _, s = fireball(N, 20.0, 4.0)
        push!(firsts, abs(s[end][2, 1])/maximum(abs, s[end][2, :]))
    end
    shrinks = all(diff(firsts) .< 1e-12)
    gate("A2", firsts[end] < 0.05 && shrinks,
        @sprintf("u^r → 0 at the axis: |u^r|(first cell)/max|u^r| = %s over N = %s, monotone",
                 join((@sprintf("%.2e", x) for x in firsts), " → "), join(Ns, "/")))

    # Self-convergence of T: a first-order upwind scheme, so the deviation between successive grids
    # should fall roughly like 1/N. Reported as an order, gated loosely — the number that matters is
    # that it CONVERGES, not its exact value at these resolutions.
    # One ghost cell per side; cell centres sit at (i-½)Δr, so a 2N grid's cells 2i-1 and 2i tile
    # the N grid's cell i EXACTLY. Block-average the fine grid down rather than sub-sampling it —
    # sub-sampling compares different points and reports the profile's slope as an error.
    interior(a) = a[2:end-1]
    sols = [interior(fireball(N, 20.0, 3.0)[2][end][1, :]) for N in (100, 200, 400)]
    coarsen(a) = [(a[2i-1] + a[2i])/2 for i in 1:length(a)÷2]
    l2(a, b)   = sqrt(sum((a .- b).^2)/length(a))
    d12 = l2(sols[1], coarsen(sols[2])); d23 = l2(sols[2], coarsen(sols[3]))
    p   = log2(d12/max(d23, 1e-300))
    gate("A3", d23 < d12 && p > 0.5,
        @sprintf("T self-converges on the fireball: ‖Δ‖ %.3e → %.3e over N = 100/200/400, observed order %.2f", d12, d23, p))
end
gate_A()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("C · cost")
# ════════════════════════════════════════════════════════════════════════════════════════════════
function gate_C()
    ts = Float64[]
    for N in (100, 200, 400)
        fireball(N, 20.0, 3.0)                       # warm-up: the first call carries compilation
        t = @elapsed fireball(N, 20.0, 3.0)
        push!(ts, t)
        @printf("        N_r = %4d   %6.2f s   (τ 0.4 → 3)\n", N, t)
    end
    @printf("        scaling N = 100 → 400: %.2fx wall time for 4x the cells\n", ts[end]/ts[1])
    note("Deliberately NOT normalised per cell: the STEP COUNT is not fixed across these grids, so")
    note("a µs/cell figure divides by the wrong thing and reads as super-linear cost that is not there.")
    note("The kernel is EoS-bound (FLUIDUM_FULL_AUDIT.md §4) — the `thermodynamic` lookups inside the")
    note("matrix build dominate, not the linear algebra around them.")
end
gate_C()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("summary")
@printf("  %d/%d gates passed\n", L.p, L.p + L.f)
L.f == 0 || println("  FAILED: ", join(L.names, ", "))
exit(L.f == 0 ? 0 : 1)
