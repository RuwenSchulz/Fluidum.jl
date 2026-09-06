#!/usr/bin/env julia
# ════════════════════════════════════════════════════════════════════════════════════════════════
# bench/bench_2p1d_viscous.jl — THE 2+1D VISCOUS SOLVER AGAINST CLOSED FORMS.
#
# `matrix2d_visc!` (src/Matrix/2d_viscous.jl) evolves
#
#       φ = (T, u^x, u^y, π^{yy}, π^{zz}, π^{xy}, Π)
#
# on a transverse-Cartesian, boost-invariant Milne slice.  Until 2026-09-03 its entire test
# coverage was one `@test all(isfinite.(result))` in `test/runtests.jl`: the solver had never been
# compared against a closed form, a limit, or another solver.  This file is that comparison, at the
# level of the PDE — every reference here is a closed form or an ODE integrated inside this file.
#
#   M  characteristic structure, from the matrix alone (closed-form speeds, causal envelope)
#   B  BJORKEN     a uniform state must integrate the 0+1D Israel–Stewart equations EXACTLY
#   S  SOUND       phase velocity and attenuation against the exact MIS dispersion relation
#   G  GUBSER      the ideal conformal solution, evolved in 2-D
#   A  STRUCTURE   axisymmetry, and the 2-D solve against the 1-D cylindrical solver
#   C  COST        cells/s and the scaling
#
#   julia --project=Julia Julia/Fluidum.jl/bench/bench_2p1d_viscous.jl
#   FLUIDUM_2DV_FULL=1 …   the finer grids in S, G and A (~20 min instead of ~6)
#
# COMPANION FILE.  `Julia/Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl` (a different
# session, same day) carries the DERIVATION side: `Julia/tools/derive_2p1d_viscous.wls` regenerates
# both the 1-D cylindrical and the 2-D Cartesian matrices from scratch, reproduces the shipped 1-D
# production kernel entry by entry (7.8e-16), and audits the shipped 2-D matrix against the derived
# one.  Read the two together: that file says which ENTRIES are wrong, this one says what the
# SOLVER does about it.
#
# 🔴 WHAT THE TWO FILES AGREE ON — the π^{xy} row is driven by the ONE-SIDED velocity gradient.
#    The shipped row carries −2η ∂_x u^y where Israel–Stewart asks for −2ησ^{xy} =
#    −η(∂_x u^y + ∂_y u^x); `A_y[π^{xy}, u^x]` is identically zero at every u and every π, and
#    d(A_x[π^{xy},u^y])/dη = 2(1+(u^x)²) is exactly twice the correct coefficient.  The audit finds
#    the shipped and derived matrices differ in exactly TEN entries, all in the two shear rows'
#    velocity columns; rows 1, 2, 3 (conservation), row 5 (π^η_η), row 7 (Π) and the whole source
#    vector are correct to round-off.
#    * It is EXACT for irrotational transverse flow (∂_x u^y = ∂_y u^x) — which is why every smooth
#      fireball started from rest has looked right, the FiVoHydro comparison included.
#    * It is wrong by −2η ω_{xy} wherever the transverse flow carries vorticity.
#    * It makes the characteristic structure ANISOTROPIC (gate M2 below).
#    Nothing in B, S, G or A is affected: Bjorken and sound-along-x and the ideal sector never touch
#    the offending entries, and an axisymmetric flow is irrotational.  Those gates therefore hold
#    for the shipped matrix and are the reference the corrected matrix must also meet.
#
# CONVENTIONS, pinned by B1 and used throughout
# ─────────────────────────────────────────────
#   π^{zz} is the MIXED Milne component π^η_η = τ²π^{ηη}; π^{yy} = π^y_y.  Same mixed convention as
#   `viscous_1d`'s `piphiphi`/`pietaeta` (checked in A0).  The relaxation is PLAIN ISRAEL–STEWART:
#   no δ_ππ, no τ_ππ, no λ_πΠ.  Navier–Stokes fixed points on a Bjorken background:
#       π^η_η → −4η/(3τ),   π^y_y → +2η/(3τ),   Π → −ζ/τ.
#   B1 pins that to 1e-13 on a state carrying |π^η_η| ≳ e/2, which is the sharpest statement
#   available: a δ_ππ of the size FiVoHydro carries (δ_ππ = (4/3)τ_π) would show at the 10 % level.
# ════════════════════════════════════════════════════════════════════════════════════════════════
using Fluidum, LinearAlgebra, Printf, OrdinaryDiffEq
const F    = Fluidum
const FULL = get(ENV, "FLUIDUM_2DV_FULL", "0") == "1"

mutable struct Led; p::Int; f::Int; names::Vector{String}; end
const L = Led(0, 0, String[])
function gate(id, ok::Bool, msg)
    ok ? (L.p += 1) : (L.f += 1; push!(L.names, id))
    @printf("  [%s] %-5s %s\n", ok ? "PASS" : "FAIL", id, msg); ok
end
note(m)    = println("        ", m)
section(t) = (println(); println("── ", t, " ", "─"^max(0, 96 - length(t))))

const EOS   = FluiduMEoS()
const CONF  = F.ConformalEOS(40.0)
const SHEAR = QGPViscosity(0.1, 0.2)                 # η/s = 0.1, τ_π = η/(s T C_s) = 5(η/s)/T̂
const BULK  = SimpleBulkViscosity(0.083, 15.0)
const PAR   = F.FluidProperties(EOS,  SHEAR,          BULK,                ZeroDiffusion())
const PAR_I = F.FluidProperties(EOS,  ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())
const PARCI = F.FluidProperties(CONF, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())

prs(T, eos = EOS)  = F.pressure(T, eos)
ent(T, eos = EOS)  = F.pressure_derivative(T, Val(1), eos)
dst(T, eos = EOS)  = F.pressure_derivative(T, Val(2), eos)
ene(T, eos = EOS)  = T*ent(T, eos) - prs(T, eos)
enth(T, eos = EOS) = T*ent(T, eos)
cs2(T, eos = EOS)  = ent(T, eos)/(T*dst(T, eos))
"(η, τ_π, ζ, τ_Π) exactly as `matrix2d_visc!` builds them"
function coeffs(T, par = PAR)
    th = thermodynamic(T, par.eos)
    (F.viscosity(T, th, par.shear), F.τ_shear(T, th, par.shear),
     F.bulk_viscosity(T, th, par.bulk), F.τ_bulk(T, th, par.bulk))
end
"(A_t⁻¹A_x, A_t⁻¹A_y, A_t⁻¹s) at one state"
function sys(T, ux, uy, pyy, pzz, pxy, PiB, τ, par = PAR)
    A = [zeros(7,7), zeros(7,7)]; S = zeros(7)
    F.matrix2d_visc!(A, S, [T,ux,uy,pyy,pzz,pxy,PiB], τ, (1.0,1.0), par)
    A[1], A[2], S
end
spd(Ax, Ay, θ) = eigvals(cos(θ)*Ax + sin(θ)*Ay)
grid2d(N, Lx, M = N, Ly = Lx) = DiscreteFields(F.viscous_2d(),
    CartesianDiscretization(F.SymmetricInterval(N, Lx), F.SymmetricInterval(M, Ly)), Float64)
cells(N, Lx) = [-Lx + (i-0.5)*2Lx/N for i in 1:N]

println("="^100)
println("Fluidum 2+1D VISCOUS HYDRO — closed-form validation ladder   (", FULL ? "FULL" : "fast", " mode)")
println("="^100)

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("M · the characteristic structure — closed forms, no PDE solve")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# The solver stores THREE shear components; the other four are fixed by u_μπ^{μν} = 0 and π^μ_μ = 0:
#   π^{xx} = [ −(π^{yy}+π^{zz})(u^τ)² + 2u^xu^yπ^{xy} + (u^y)²π^{yy} ] / (1 + (u^y)²)
#   π^{τν} = (u^xπ^{xν} + u^yπ^{yν})/u^τ
# `1 + (u^y)²` is the denominator that runs through the whole generated matrix — M1 checks that this
# really is the closure, by verifying the reconstructed tensor is traceless and u-orthogonal.
function pi_full(ux, uy, pyy, pzz, pxy)
    ut  = sqrt(1 + ux^2 + uy^2)
    pxx = (-(pyy + pzz)*ut^2 + 2ux*uy*pxy + uy^2*pyy)/(1 + uy^2)
    ptx = (ux*pxx + uy*pxy)/ut
    pty = (ux*pxy + uy*pyy)/ut
    ptt = (ux*ptx + uy*pty)/ut
    # (τ, x, y, z) with the local Cartesian z = τη at midrapidity, so g = diag(−1, 1, 1, 1)
    [ptt ptx pty 0.0; ptx pxx pxy 0.0; pty pxy pyy 0.0; 0.0 0.0 0.0 pzz]
end
function gate_M1()
    g = Diagonal([-1.0, 1.0, 1.0, 1.0]); wtr = 0.0; wor = 0.0
    for ux in (-2.0,-0.5,0.0,0.7,3.0), uy in (-1.5,0.0,0.3,2.0),
        (pyy,pzz,pxy) in ((0.1,-0.2,0.05), (-0.03,0.07,-0.11), (0.5,0.5,0.5))
        Π = pi_full(ux, uy, pyy, pzz, pxy); u = [sqrt(1+ux^2+uy^2), ux, uy, 0.0]
        sc = max(maximum(abs, Π), 1e-30)
        wtr = max(wtr, abs(tr(g*Π))/sc); wor = max(wor, maximum(abs, (u'*g*Π)')/sc)
    end
    gate("M1", wtr < 1e-13 && wor < 1e-13,
        @sprintf("the 3 stored components close onto a traceless, u-orthogonal π^{μν} (|π^μ_μ| %.1e, |u_μπ^{μν}| %.1e, relative)", wtr, wor))
end
gate_M1()

# M2 — THE CANARY.  Characteristic speeds are properties of the system, not of the variables used
# to write it: rotate the state and the probe direction together and the seven speeds must not move.
#
# Run TWICE, once against each matrix, because as of 2026-09-03 `matrix2d_visc!` can be routed to
# either: the shipped generated matrix (default) or the one `Julia/tools/derive_2p1d_viscous.wls`
# derives from scratch (`FLUIDUM_2D_DERIVED=1`).  The gate asserts on the DERIVED matrix — that is
# the statement about the physics — and reports the shipped one alongside, so the size of the defect
# and the fact that the flag really cures it are both visible in one place.
function isotropy_residual()
    worst_ideal = 0.0; worst_all = 0.0
    for T in (0.15,0.30,0.50), U in (0.0,0.6,1.3,2.5), fr in (0.0, 0.12), τ in (0.8,4.0)
        e = ene(T); pzz = -fr*e; pyy = -pzz/2; pxy = 0.3fr*e; PiB = -0.3fr*e
        Ax0, Ay0, _ = sys(T, U, 0.0, pyy, pzz, pxy, PiB, τ)
        for ψ in range(0, 2π, length=9)[1:8]
            c, s = cos(ψ), sin(ψ)
            Pin  = [c -s; s c]*pi_full(U,0.0,pyy,pzz,pxy)[2:3,2:3]*[c s; -s c]
            Ax, Ay, _ = sys(T, U*c, U*s, Pin[2,2], pzz, Pin[1,2], PiB, τ)
            for θ in range(0, 2π, length=7)[1:6]
                d = maximum(abs, sort(real.(spd(Ax,Ay,θ+ψ))) .- sort(real.(spd(Ax0,Ay0,θ))))
                worst_all = max(worst_all, d); fr == 0.0 && (worst_ideal = max(worst_ideal, d))
            end
        end
    end
    (worst_all, worst_ideal)
end
function gate_M2()
    had = F.VISC_2D_DERIVED[]
    F.VISC_2D_DERIVED[] = false; sh_all, sh_id = isotropy_residual()
    F.VISC_2D_DERIVED[] = true;  dv_all, dv_id = isotropy_residual()
    F.VISC_2D_DERIVED[] = had
    gate("M2", dv_all < 1e-11,
        @sprintf("rotational isotropy of the 7 characteristic speeds, DERIVED matrix: %.2e (%.2e with π = Π = 0)", dv_all, dv_id))
    note(@sprintf("the SHIPPED matrix gives %.3f c (%.3f c with π = Π = 0) — the π^{xy} defect in the header.", sh_all, sh_id))
    note("on the SHIPPED matrix the transverse shear channel propagates at √(2C_s) along x and at 0")
    note("along y — a factor √2 too fast, see M3b; the ideal sector is unaffected, which is why B,")
    note("G and A below hold for either matrix.")
    note("FLUIDUM_2D_DERIVED selects the corrected matrix everywhere, including the 10-field")
    note("`matrix2d_visc_HQ_BG!`.  DEFAULT ON since 2026-09-03 (`a1d9f924`); set it to 0 only to")
    note("reproduce results already on disk.  It is not slower (measured 170 ns vs 193 ns).")
end
gate_M2()

function gate_M3()
    worst = 0.0
    for T in (0.15,0.20,0.30,0.45,0.60)
        η, τπ, ζ, τΠ = coeffs(T); w = enth(T)
        vf = sqrt(cs2(T) + (4/3)*η/(w*τπ) + ζ/(w*τΠ))
        vm = maximum(abs, real.(eigvals(sys(T,0.0,0.0,0.0,0.0,0.0,0.0,2.0)[1])))
        worst = max(worst, abs(vm - vf)/vf)
        note(@sprintf("T = %.2f GeV   v_front = %.8f   closed form %.8f", T, vm, vf))
    end
    gate("M3", worst < 1e-12,
        @sprintf("at rest the largest characteristic ≡ √(c_s² + (4/3)η/(wτ_π) + ζ/(wτ_Π)), the MIS front (%.2e)", worst))
    # 🔴 CORRECTED 2026-09-03 (evening), BY MEASUREMENT — this gate used to assert √(2C_s), and
    # that number was the DEFECT, not the physics.  It was written while the shipped matrix was the
    # default, and the shipped π^{xy} row carries −2η ∂_x u^y where Israel–Stewart asks for
    # −2ησ^{xy} = −η(∂_x u^y + ∂_y u^x) (the header's "exactly twice the correct coefficient").  For
    # a perturbation depending on x alone that doubles the u^y ↔ π^{xy} coupling, hence doubles the
    # squared speed — so the legacy matrix propagates the transverse shear channel at √(2C_s) and
    # the corrected one at √(C_s), which is the Israel–Stewart shear-channel velocity
    # v² = η/(τ_π(e+P)).  Measured on the derived matrix: the eigenvalue is √(C_s) EXACTLY
    # (0.316228 at C_s = 0.1, 0.447214 at 0.2, 0.591608 at 0.35, both temperatures, residual 3e-16
    # once the reference is right), and on the legacy matrix it is √(2C_s) to 3.3e-16.
    #
    # This is the failure class B19 already recorded — a gate asserting a hand copy of the formula
    # under test, and so certifying the bug.  It survived the 2026-09-03 audit because the flag was
    # still default OFF when the audit ran; the default flip (`a1d9f924`) is what exposed it, and
    # `run_suite.jl bench` is what ran the file afterwards.  M3 above is untouched and still passes
    # at 8e-16, which is the independent check that the shear NORMALISATION is otherwise right: the
    # sound front carries (4/3)η/(wτ_π) through a different row.
    w2 = 0.0
    for Cs in (0.1,0.2,0.35), T in (0.2,0.45)
        pr = F.FluidProperties(EOS, QGPViscosity(0.1,Cs), ZeroBulkViscosity(), ZeroDiffusion())
        ev = sort(real.(eigvals(sys(T,0.0,0.0,0.0,0.0,0.0,0.0,2.0,pr)[1])))
        want = F.VISC_2D_DERIVED[] ? sqrt(Cs) : sqrt(2Cs)   # the legacy matrix's speed is the defect
        w2 = max(w2, minimum(abs.(ev .- want)))
    end
    gate("M3b", w2 < 1e-12,
        @sprintf("the transverse shear channel sits at √(η/(wτ_π)) = √(C_s), independent of T and of the EoS (%.2e)%s",
                 w2, F.VISC_2D_DERIVED[] ? "" : "  [LEGACY matrix: asserted against √(2C_s), which is the defect]"))
end
gate_M3()

function gate_M4()
    worst = 0.0
    for T in (0.15,0.30,0.50), U in (0.0,0.4,1.2,3.0), τ in (1.0,5.0)
        c = sqrt(cs2(T)); v = U/sqrt(1+U^2)
        ev = sort(real.(eigvals(sys(T,U,0.0,0.0,0.0,0.0,0.0,τ,PAR_I)[1])))
        worst = max(worst, abs(ev[end] - (v+c)/(1+v*c))/abs((v+c)/(1+v*c)),
                           abs(ev[1]   - (v-c)/(1-v*c))/max(abs((v-c)/(1-v*c)), 1e-3))
    end
    gate("M4", worst < 1e-11,
        @sprintf("with viscosity off the extreme characteristics ≡ the sound cone (v±c_s)/(1±v c_s) (%.2e)", worst))
end
gate_M4()

# M5 — the causal envelope.  MIS is hyperbolic and subluminal only inside a bounded region of
# (π, Π) and Fluidum carries no regulator, so where that boundary lies IS the operating envelope.
# Measured, not asserted; the gate only asks that it contain what a fireball reaches.
function envelope(T, U, τ; mode = :both, fmax = 1.0, df = 0.005)
    e = ene(T)
    for fr in 0.0:df:fmax
        pzz = mode === :bulk  ? 0.0 : -fr*e
        PiB = mode === :shear ? 0.0 : -fr*e
        Ax, Ay, _ = sys(T, U, 0.0, -pzz/2, pzz, 0.0, PiB, τ)
        for θ in range(0,2π,length=17)[1:16]
            ev = spd(Ax, Ay, θ)
            (maximum(abs, imag.(ev)) > 1e-6 || maximum(abs, real.(ev)) > 1.0) && return fr
        end
    end
    fmax
end
function gate_M5()
    note("the dissipative fraction at which the system first turns acausal or non-hyperbolic:")
    note(string(rpad("T",7), rpad("|u|",6), rpad("shear only",13), rpad("bulk only",13), "both"))
    worst = Inf
    for T in (0.15,0.20,0.30,0.45), U in (0.0,1.0,2.0)
        a = envelope(T,U,2.0; mode=:shear); b = envelope(T,U,2.0; mode=:bulk); c = envelope(T,U,2.0)
        note(@sprintf("%-7.2f%-6.1f%-13.3f%-13.3f%.3f", T, U, a, b, c)); worst = min(worst, c)
    end
    gate("M5", worst > 0.10,
        @sprintf("the causal envelope is at worst |π|/e = %.3f, wider than the ≲0.06 a fireball reaches at these T", worst))
    note("bulk is the binding constraint, and it binds NEAR THE LORENTZIAN PEAK OF ζ(T), not at low T:")
    note("the threshold tracks ζ/(w τ_Π) inversely and that ratio is peaked at T = 0.175, so the")
    note("envelope is NARROWEST there (|Π|/e ≈ 0.41) and widens on BOTH sides — ≈0.635 at T = 0.03,")
    note("≈0.795 at T = 0.60. ⚠ Reading it as a cold-tail constraint is a trap: it cost this session")
    note("a wrong diagnosis of a coarse-grid artefact (examples/2p1d_viscous/03_elliptic_fireball.jl).")
    note("Shear alone stays subluminal to |π|/e ≥ 1.0 at rest at every T tested.")
end
gate_M5()

function gate_M6()
    bad = 0.0; Tb = 0.0
    for T in (0.02,0.05,0.10,0.15,0.30,0.60)
        th = thermodynamic(T, EOS); η, τπ, ζ, τΠ = coeffs(T)
        At = F.two_d_viscous_matrix([T,0.3,-0.2,0.01,-0.02,0.003,0.005], 2.0, th.pressure,
                th.pressure_derivative[1], th.pressure_hessian[1], ζ, η, τπ, τΠ)[1]
        κ = cond(Matrix(At)); κ > bad && (bad = κ; Tb = T)
        note(@sprintf("T = %.2f GeV   cond(A_t) = %.3e", T, κ))
    end
    gate("M6", bad < 1e12, @sprintf("A_t stays invertible down to T = 0.02 GeV (worst cond %.2e at T = %.2f)", bad, Tb))
end
gate_M6()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("B · Bjorken — the nonlinear viscous sector, exactly")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# A transversely uniform state has no transverse gradient, so the 2-D solver must integrate the
# 0+1D Israel–Stewart system.  That system is written here from the physics, not from the matrix:
#   De = −(e + P + Π + π^η_η)/τ ,  τ_π Dπ^η_η = −π^η_η − 4η/(3τ) ,  τ_Π DΠ = −Π − ζ/τ
# and de/dT = T d²P/dT².  This is the ONLY gate that tests the shear sector nonlinearly against a
# reference written independently of the solver, and it is what pins the conventions.
function bjorken_ref(T0, τ0, τ1; par = PAR)
    function rhs!(du, U, p, τ)
        T, pe, PB = U
        η, τπ, ζ, τΠ = coeffs(T, par)
        du[1] = (-(ene(T) + prs(T) + PB + pe)/τ)/(T*dst(T))
        du[2] = (-pe - 4η/(3τ))/τπ
        du[3] = (-PB - ζ/τ)/τΠ
    end
    solve(ODEProblem(rhs!, [T0,0.0,0.0], (τ0,τ1)), Tsit5(); reltol=1e-12, abstol=1e-14)(τ1)
end
function gate_B()
    τ0, τ1, T0 = 0.4, 3.0, 0.35
    d   = grid2d(24, 12.0)
    phi = set_array((x,y)->T0, :temperature, d)
    for f in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    u   = F.oneshoot(d, F.matrix2d_visc!, PAR, phi, (τ0,τ1); reltol=1e-11, abstol=1e-13)(τ1)
    r   = bjorken_ref(T0, τ0, τ1)
    ic, jc = 13, 13
    eT  = abs(u[1,ic,jc] - r[1])/abs(r[1])
    ez  = abs(u[5,ic,jc] - r[2])/abs(r[2])
    eB  = abs(u[7,ic,jc] - r[3])/abs(r[3])
    gate("B1", max(eT,ez,eB) < 1e-10,
        @sprintf("a uniform state integrates the 0+1D MIS equations: T %.2e, π^η_η %.2e, Π %.2e (relative)", eT, ez, eB))
    note(@sprintf("at τ = %.1f: T = %.8f, π^η_η = %.8f (= %.2f×e), Π = %.8f", τ1, r[1], r[2], r[2]/ene(r[1]), r[3]))
    note(@sprintf("a δ_ππ = (4/3)τ_π term of FiVoHydro's size would move π^η_η by ≈ %.1f %% on this state — it is absent",
                  100*abs((4/3)*coeffs(r[1])[2]/τ1)))
    # B2 — tracelessness is DYNAMICAL here: π^{yy} and π^{zz} are evolved by separate rows
    tr = abs(u[4,ic,jc] + u[5,ic,jc]/2)/abs(u[5,ic,jc])
    gate("B2", tr < 1e-12,
        @sprintf("π^y_y = −π^η_η/2 emerges from two independently evolved rows (residual %.2e)", tr))
    # B3 — transverse uniformity and rest must be preserved exactly
    sp = maximum(abs, u[1,2:25,2:25] .- u[1,ic,jc])/abs(u[1,ic,jc])
    mv = maximum(abs, u[2,2:25,2:25]) + maximum(abs, u[3,2:25,2:25])
    gate("B3", sp < 1e-14 && mv < 1e-14,
        @sprintf("uniformity and rest are preserved to round-off (spread %.1e, max|u| %.1e)", sp, mv))
    # B4 — ideal Bjorken: τ·s is conserved
    phi = set_array((x,y)->T0, :temperature, d)
    for f in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    ui  = F.oneshoot(d, F.matrix2d_visc!, PAR_I, phi, (τ0,τ1); reltol=1e-11, abstol=1e-13)(τ1)
    es  = abs(τ1*ent(ui[1,ic,jc]) - τ0*ent(T0))/(τ0*ent(T0))
    gate("B4", es < 1e-10, @sprintf("with viscosity off, τ·s is conserved on the Bjorken solution (%.2e)", es))
end
gate_B()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("S · sound — the viscous sector against the exact MIS dispersion relation")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# Linearising the Israel–Stewart system about a static uniform state (written from the physics, not
# from the matrix) gives, for a plane wave along x,
#
#     ω² = k²c_s² − i ω k² Z(ω)/w ,     Z(ω) = (4/3)η/(1 − iωτ_π) + ζ/(1 − iωτ_Π) ,   w = e + P
#
# whose small-k limit is ω ≈ c_s k − (i/2)k²((4/3)η + ζ)/w.  The gate solves the FULL root, because
# at the k used here ωτ_π ≈ 0.3 and the Navier–Stokes limit is already 10 % off.
#
# ⚠ δ_ππ π θ is SECOND order in the perturbation (π and θ are both first order), so this gate is
# blind to δ_ππ by construction — which is exactly why it is a clean test of the first-order sector.
# The run sits at τ₀ = 200 fm/c so the Milne source terms (∝1/τ) are ~0.5 % of ω, and the amplitude
# is extracted from the DIFFERENCE of a perturbed and an unperturbed run so the residual background
# drift cancels.  The seed is the exact right-moving eigenvector, so the signal decays monotonically
# instead of beating against a counter-propagating twin.
function mis_omega(k, T, par = PAR)
    η, τπ, ζ, τΠ = coeffs(T, par); w = enth(T); c2 = cs2(T)
    Z(ω) = (4/3)*η/(1 - im*ω*τπ) + ζ/(1 - im*ω*τΠ)
    f(ω) = ω^2 - k^2*c2 + im*ω*k^2*Z(ω)/w
    ω = sqrt(c2)*k - 0.5im*k^2*((4/3)*η + ζ)/w
    for _ in 1:200
        h = 1e-7*max(abs(ω), 1e-3); ω -= f(ω)/((f(ω+h) - f(ω-h))/(2h))
    end
    ω
end
function sound_eigvec(k, T, par = PAR)
    ω = mis_omega(k, T, par); η, τπ, ζ, τΠ = coeffs(T, par); w = enth(T)
    v   = ω/(k*w)                                   # normalised to δe = 1
    pxx = -(4/3)*η*(im*k*v)/(1 - im*ω*τπ)
    (ω, 1/(T*dst(T)), v, -pxx/2, -pxx/2, -ζ*(im*k*v)/(1 - im*ω*τΠ))
end
function sound_measure(T0, k, Lx, N, τ0, τ1; xmax = 12.0, nout = 17, amp = 1e-4, par = PAR)
    ω, dT, v, pyy, pzz, PiB = sound_eigvec(k, T0, par); s = amp*T0/abs(dT)
    d = grid2d(N, Lx, 4, 2.0); τs = collect(range(τ0, τ1; length=nout)); res = Dict()
    for pert in (false, true)
        q = pert ? s : 0.0
        phi = set_array((x,y)->T0 + q*real(dT*exp(im*k*x)), :temperature, d)
        set_array!(phi, (x,y)->q*real(v  *exp(im*k*x)), :ux,   d)
        set_array!(phi, (x,y)->q*real(pyy*exp(im*k*x)), :piyy, d)
        set_array!(phi, (x,y)->q*real(pzz*exp(im*k*x)), :pizz, d)
        set_array!(phi, (x,y)->q*real(PiB*exp(im*k*x)), :piB,  d)
        set_array!(phi, (x,y)->0.0, :uy, d); set_array!(phi, (x,y)->0.0, :pixy, d)
        sol = F.oneshoot(d, F.matrix2d_visc!, par, phi, (τ0,τ1); reltol=1e-11, abstol=1e-13)
        for τ in τs; res[(pert,τ)] = sol(τ); end
    end
    xs = cells(N, Lx); sel = findall(x->abs(x) <= xmax, xs)
    amps = ComplexF64[]
    for τ in τs
        δ = [res[(true,τ)][1,i+1,3] - res[(false,τ)][1,i+1,3] for i in 1:N]
        push!(amps, 2*sum(δ[i]*exp(-im*k*xs[i]) for i in sel)/length(sel))
    end
    lm = log.(abs.(amps)); ph = angle.(amps)
    for i in 2:length(ph)
        while ph[i]-ph[i-1] >  π; ph[i] -= 2π; end
        while ph[i]-ph[i-1] < -π; ph[i] += 2π; end
    end
    X = hcat(ones(nout), τs .- τ0)
    (-(X\ph)[2], -(X\lm)[2], ω)
end
function gate_S()
    T0 = 0.30; Ns = FULL ? (100,200,400,800) : (100,200,400)
    note("λ = 6 fm (k = 1.047 fm⁻¹), L = 20 fm, τ = 200 → 208 fm/c, amplitude 1e-4:")
    note(string(rpad("N",7), rpad("Δx",9), rpad("Re ω",12), rpad("rel",10), rpad("Im ω",12), rpad("excess",11), "order"))
    k = 2π/6.0; er = Float64[]; ei = Float64[]; last = (0.0,0.0,0.0+0im)
    for N in Ns
        wr, wi, ω = sound_measure(T0, k, 20.0, N, 200.0, 208.0)
        push!(er, abs(wr - real(ω))); push!(ei, wi + imag(ω)); last = (wr, wi, ω)
        note(@sprintf("%-7d%-9.4f%-12.6f%-10.2e%-12.6f%-11.6f%s", N, 40.0/N, wr, er[end]/abs(real(ω)), wi, ei[end],
                      length(ei) < 2 ? "—" : @sprintf("%.2f / %.2f", er[end-1]/er[end], ei[end-1]/ei[end])))
    end
    ω = last[3]
    gate("S1", er[end]/abs(real(ω)) < 1e-2 && er[1] > er[end],
        @sprintf("the phase velocity converges onto the exact MIS root: Re ω = %.6f vs %.6f (%.2e), from %.2e at the coarsest grid",
                 last[1], real(ω), er[end]/abs(real(ω)), er[1]/abs(real(ω))))
    ords = [ei[i]/ei[i+1] for i in 1:length(ei)-1]
    gate("S2", all(o -> 1.7 < o < 2.3, ords),
        @sprintf("the attenuation excess is FIRST ORDER in Δx — %s per doubling — i.e. it is the upwind scheme's numerical viscosity, not a missing term",
                 join([@sprintf("%.2f", o) for o in ords], ", ")))
    Dn = 2*ei[end]/k^2
    η, τπ, ζ, τΠ = coeffs(T0); Γs = ((4/3)*η + ζ)/enth(T0)
    rich = -imag(ω) + 2ei[end] - ei[end-1]
    note(@sprintf("Richardson from the two finest grids: Im ω → %.6f against the exact %.6f (%.1f %%)",
                  rich, -imag(ω), 100*abs(rich + imag(ω))/abs(imag(ω))))
    Dprev = 2*ei[end-1]/k^2
    gate("S3", abs(Dn/(40.0/Ns[end]) - Dprev/(80.0/Ns[end]))/(Dn/(40.0/Ns[end])) < 0.10,
        @sprintf("NUMERICAL VISCOSITY: D_num = %.4f fm at Δx = %.3f fm, i.e. %.2f·Δx — and %.2f·Δx one grid coarser, so the coefficient is a property of the scheme, not of the run",
                 Dn, 40.0/Ns[end], Dn/(40.0/Ns[end]), Dprev/(80.0/Ns[end])))
    note(@sprintf("the physical Γ_s = ((4/3)η + ζ)/(e+P) = %.4f fm at T = %.2f GeV, so the scheme's own", Γs, T0))
    note(@sprintf("viscosity equals the physical one at Δx ≈ %.3f fm and DOMINATES above it.", Γs/(Dn/(40.0/Ns[end]))))
    note("⚠ the FiVoHydro comparison (COMPARISON_2P1D.md) ran its finest grid at Δx = 0.15 fm.")
end
gate_S()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("G · Gubser — the ideal conformal solution, evolved in 2-D")
# ════════════════════════════════════════════════════════════════════════════════════════════════
const QG = 1/4.3
gT(τ, r; T0=1.2) = T0*(2QG)^(2/3)/(τ^(1/3)*(1 + 2QG^2*(τ^2+r^2) + QG^4*(τ^2-r^2)^2)^(1/3))
gu(τ, r) = (v = 2QG^2*τ*r/(1 + QG^2*(τ^2+r^2)); v/sqrt(1-v^2))
function gubser_run(N, Lx, τ0, τ1; rcut = 6.0)
    d   = grid2d(N, Lx)
    phi = set_array((x,y)->gT(τ0, hypot(x,y)), :temperature, d)
    set_array!(phi, (x,y)->(r = hypot(x,y); r < 1e-12 ? 0.0 : gu(τ0,r)*x/r), :ux, d)
    set_array!(phi, (x,y)->(r = hypot(x,y); r < 1e-12 ? 0.0 : gu(τ0,r)*y/r), :uy, d)
    for f in (:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    u  = F.oneshoot(d, F.matrix2d_visc!, PARCI, phi, (τ0,τ1); reltol=1e-9, abstol=1e-11)(τ1)
    xs = cells(N, Lx); eT = 0.0; eU = 0.0; asym = 0.0
    for i in 1:N, j in 1:N
        r = hypot(xs[i], xs[j]); r > rcut && continue
        eT = max(eT, abs(u[1,i+1,j+1] - gT(τ1,r))/gT(τ1,r))
        ur = gu(τ1,r); eU = max(eU, abs(u[2,i+1,j+1] - (r<1e-9 ? 0.0 : ur*xs[i]/r))/max(abs(ur),1e-3))
        asym = max(asym, abs(u[1,i+1,j+1] - u[1,j+1,i+1]))
    end
    (eT, eU, asym)
end
function gate_G()
    cf = 0.0
    for T in (0.2,0.5,1.0,1.5); cf = max(cf, abs(ene(T,CONF)/(3prs(T,CONF)) - 1)); end
    gate("G0", cf < 1e-12, @sprintf("the conformal EoS Gubser needs really is conformal: max |e/3P − 1| = %.2e", cf))
    Ns = FULL ? (100,200,400) : (100,200)
    note("ideal Gubser (q = 1/4.3 fm⁻¹), τ = 1 → 2 fm/c, box ±12 fm, compared over r ≤ 6 fm:")
    note(string(rpad("N",7), rpad("Δx",9), rpad("max rel dev T",16), rpad("max rel dev u^x",17), rpad("x↔y",11), "order"))
    prev = NaN; ok = true; last = 0.0; ords = Float64[]
    for N in Ns
        eT, eU, a = gubser_run(N, 12.0, 1.0, 2.0)
        isnan(prev) || push!(ords, log2(prev/eT))
        note(@sprintf("%-7d%-9.4f%-16.4e%-17.4e%-11.2e%s", N, 24.0/N, eT, eU, a,
                      isnan(prev) ? "—" : @sprintf("%.2f", log2(prev/eT))))
        a > 1e-12 && (ok = false); prev = eT; last = eT
    end
    gate("G1", ok && last < 5e-3 && all(o -> o > 0.85, ords),
        @sprintf("the 2-D solver reproduces ideal Gubser and converges at first order (%s), x↔y symmetric to round-off",
                 join([@sprintf("%.2f", o) for o in ords], ", ")))
    note("first order is the expected rate: the spatial scheme is first-order upwind (see S2).")
    note("the ideal sector is x↔y symmetric at round-off — the π^{xy} defect lives only in the shear rows.")
end
gate_G()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("A · structure — the 2-D viscous solve against the 1-D cylindrical solver")
# ════════════════════════════════════════════════════════════════════════════════════════════════
# `matrix1d_visc!` is the independently validated 1-D cylindrical kernel (gates S10/S11b of
# Projects/FluidumValidation/bench_fluidum_full.jl, and the matrix the derivation reproduces entry
# by entry).  An axisymmetric 2-D run must converge onto it.  This is the only gate that exercises
# the full nonlinear 2-D viscous system against another solver — and axisymmetric flow is
# irrotational, so it is blind to the π^{xy} defect by construction.
T_ws(r) = 0.45/(1 + exp((r - 6.0)/0.6)) + 0.05
function gate_A0()
    w = 0.0
    for T in (0.15,0.25,0.35,0.5)
        th = thermodynamic(T,EOS); dpt = ent(T); dptt = dst(T)
        a = (F.viscosity(T,dpt,PAR.shear), F.τ_shear(T,dpt,PAR.shear),
             F.bulk_viscosity(T,dpt,PAR.bulk), F.τ_bulk(T,dpt,dptt,PAR.bulk))
        b = (F.viscosity(T,th,PAR.shear),  F.τ_shear(T,th,PAR.shear),
             F.bulk_viscosity(T,th,PAR.bulk), F.τ_bulk(T,th,PAR.bulk))
        w = max(w, maximum(abs.((a .- b)./b)))
    end
    gate("A0", w < 1e-12,
        @sprintf("the scalar-entropy (1-D) and Thermodynamic (2-D) coefficient paths give the same (η, τ_π, ζ, τ_Π) (%.2e)", w))
    note("these were DISJOINT method sets until 2026-09-02 — each entry point worked with exactly")
    note("the viscosity model the other could not take.  The gate exists so that cannot come back.")
end
gate_A0()
function run1d_visc(Nr, rmax, τ0, τ1)
    d   = DiscreteFields(F.viscous_1d(), CartesianDiscretization(F.OriginInterval(Nr, rmax)), Float64)
    phi = set_array(T_ws, :temperature, d)
    for f in (:ur,:piphiphi,:pietaeta,:piB); set_array!(phi, r->0.0, f, d); end
    (F.oneshoot(d, F.matrix1d_visc!, PAR, phi, (τ0,τ1); reltol=1e-10, abstol=1e-12)(τ1),
     [rmax*(i-0.5)/Nr for i in 1:Nr])
end
function run2d_visc(N, Lx, τ0, τ1)
    d   = grid2d(N, Lx)
    phi = set_array((x,y)->T_ws(hypot(x,y)), :temperature, d)
    for f in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    (F.oneshoot(d, F.matrix2d_visc!, PAR, phi, (τ0,τ1); reltol=1e-10, abstol=1e-12)(τ1), cells(N, Lx))
end
lerp(rs, v, r) = (i = searchsortedfirst(rs, r);
                  i <= 1 ? v[1] : i > length(rs) ? v[end] : v[i-1] + (v[i]-v[i-1])*(r-rs[i-1])/(rs[i]-rs[i-1]))
function gate_A()
    τ0, τ1, Lx = 0.6, 3.0, 16.0
    Nr = 1200; u1, rs = run1d_visc(Nr, Lx, τ0, τ1)
    note(@sprintf("1-D reference: N_r = %d (Δr = %.4f fm), Woods–Saxon T₀, at rest, τ %.1f → %.1f", Nr, Lx/Nr, τ0, τ1))
    note(string(rpad("N",7), rpad("Δx",9), rpad("T",11), rpad("u^r",11), rpad("π^y_y",11), rpad("π^η_η",11), rpad("Π",11), "azim spread"))
    Ns = FULL ? (80,160,320) : (80,160)
    prev = NaN; ords = Float64[]; lastT = 0.0; azs = Float64[]
    ref = [u1[m,1:Nr] for m in 1:5]; sc = [maximum(abs, ref[m]) for m in 1:5]
    hot  = findall(T -> T > 0.12, ref[1])
    frac = maximum(abs.(ref[4][hot])./[ene(T) for T in ref[1][hot]])
    frac_all = maximum(abs.(ref[4])./[ene(T) for T in ref[1]])
    for N in Ns
        u2, xs = run2d_visc(N, Lx, τ0, τ1); jc = N÷2
        er = zeros(5); az = 0.0
        for r in (1.0,2.0,3.0,4.0,5.0,6.0,7.0)
            i = argmin(abs.(xs .- r))
            got = [0.5*(u2[k,i+1,jc+1] + u2[k,i+1,jc+2]) for k in (1,2,4,5,7)]
            for m in 1:5; er[m] = max(er[m], abs(got[m] - lerp(rs, ref[m], r))/sc[m]); end
            vals = [u2[1, argmin(abs.(xs .- r*cos(th)))+1, argmin(abs.(xs .- r*sin(th)))+1]
                    for th in range(0,2π,length=25)[1:24]]
            az = max(az, (maximum(vals) - minimum(vals))/maximum(vals))
        end
        isnan(prev) || push!(ords, log2(prev/er[1]))
        note(@sprintf("%-7d%-9.4f%-11.3e%-11.3e%-11.3e%-11.3e%-11.3e%.3e", N, 2Lx/N, er..., az))
        prev = er[1]; lastT = er[1]; push!(azs, az)
    end
    gate("A1", lastT < 1.5e-2 && all(o -> o > 0.7, ords),
        @sprintf("the axisymmetric 2-D viscous solve converges onto the 1-D cylindrical solve in all five fields (T to %.2e at Δx = %.2f fm, order %s)",
                 lastT, 2Lx/Ns[end], join([@sprintf("%.2f", o) for o in ords], ", ")))
    note("all five fields converge together; Π is the slowest, as it is the smallest.")
    azr = [azs[i]/azs[i+1] for i in 1:length(azs)-1]
    gate("A2", azs[end] < 0.10 && all(r -> r > 1.2, azr),
        @sprintf("axisymmetry is preserved in the limit: ring spread of T = %.2e at Δx = %.2f fm, shrinking %s per doubling",
                 azs[end], 2Lx/Ns[end], join([@sprintf("%.2f", r) for r in azr], ", ")))
    note("⚠ this is the Cartesian grid's own anisotropy and it converges only at FIRST order. For")
    note("   scale, FiVoHydro's ideal-Gubser ring spread (gate G1 of its ladder) is 6.4e-4 at")
    note("   Δx = 0.05 fm and converges at SECOND order — a scheme difference, not a physics one.")
    note(@sprintf("the fireball in this gate reaches |π^η_η|/e = %.3f over T > 0.12 GeV — inside the", frac))
    note(@sprintf("   M5 envelope — but %.3f over the WHOLE grid, i.e. OUTSIDE it, in the dilute tail", frac_all))
    note("   where e → 0.  MIS is not causal there and Fluidum carries no regulator; the cells are")
    note("   cold enough that nothing is read out of them, which is a reason to mask observables at")
    note("   a temperature cut, not a reason to call the tail safe.")
end
gate_A()

# ════════════════════════════════════════════════════════════════════════════════════════════════
section("C · cost")
# ════════════════════════════════════════════════════════════════════════════════════════════════
function gate_C()
    note(string(rpad("N",7), rpad("cells",9), rpad("wall (s)", 11), "μs per cell per unit τ"))
    prev = (0, 0.0)
    for N in (80, 160)
        t = @elapsed run2d_visc(N, 16.0, 0.6, 2.0)
        note(@sprintf("%-7d%-9d%-11.2f%.3f", N, N^2, t, 1e6*t/(N^2*1.4)))
        prev = (N, t)
    end
    note("the kernel is EoS-bound and the step count is set by FLUIDUM_DTMAX, so the cost scales as N².")
    note("reported, not gated: a wall-clock assertion would fail on a loaded machine, not on a defect.")
end
gate_C()

println()
println("="^100)
@printf("%d/%d gates passed", L.p, L.p + L.f)
L.f > 0 && @printf("   FAILED: %s", join(L.names, ", "))
println(); println("="^100)
exit(L.f == 0 ? 0 : 1)
