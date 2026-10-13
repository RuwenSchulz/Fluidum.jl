# ==============================================================================================
# test/test_1p1d_viscous.jl -- the 1+1D viscous solver, `matrix1d_visc!`.
#
# THE PRODUCTION KERNEL. LangevinPaper1, AttractorHydro and SpectraDiagnostic all stand on the
# matrix asserted here, and the charm driver `matrix1d_visc_HQ!` reuses its five rows -- yet until
# 2026-09-14 it had no testset of its own. The 2+1D solver got one on 2026-09-03 (after ten wrong
# entries in its shear rows had survived a year of `@test all(isfinite.(...))`); this is its 1-D
# twin, asserting what a smoke test cannot see and staying cheap enough for every `Pkg.test()`:
# the field layout and its parities, closed-form characteristic speeds, the shear closure, the
# Bjorken limit, and the control that the shipped kernel matches an independent re-derivation.
#
# The heavy, quantitative version of the same statements (convergence on a real fireball, the
# origin-parity ladder, cost) is
#     Julia/Fluidum.jl/bench/bench_1p1d_viscous.jl
# and the audits are
#     Projects/FluidumValidation/bench_fluidum_full.jl  +  bench_fluidum_attractor.jl
# ==============================================================================================
using Test, LinearAlgebra
# no Printf: under `Pkg.test` the environment is Fluidum's deps plus [targets] test = ["Test"], and
# Printf is in neither -- importing it here would turn the suite red while still looking green
# under `--project=Julia`, whose outer environment does have it. (The trap its 2-D twin records.)
using Fluidum, OrdinaryDiffEq

_prs(T, eos)  = Fluidum.pressure(T, eos)
_ent(T, eos)  = Fluidum.pressure_derivative(T, Val(1), eos)
_dst(T, eos)  = Fluidum.pressure_derivative(T, Val(2), eos)
_ene(T, eos)  = T*_ent(T, eos) - _prs(T, eos)
_enth(T, eos) = T*_ent(T, eos)
_cs2(T, eos)  = _ent(T, eos)/(T*_dst(T, eos))

const _EOS1  = FluiduMEoS()
const _PAR1  = Fluidum.FluidProperties(_EOS1, QGPViscosity(0.1, 0.2),
                                       SimpleBulkViscosity(0.083, 15.0), ZeroDiffusion())
const _PAR1I = Fluidum.FluidProperties(_EOS1, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())

function _coeffs1(T, par = _PAR1)
    dpt  = Fluidum.pressure_derivative(T, Val(1), par.eos)
    dptt = Fluidum.pressure_derivative(T, Val(2), par.eos)
    (Fluidum.viscosity(T, dpt, par.shear),     Fluidum.τ_shear(T, dpt, par.shear),
     Fluidum.bulk_viscosity(T, dpt, par.bulk), Fluidum.τ_bulk(T, dpt, dptt, par.bulk))
end
"(A_t^-1 A_r, A_t^-1 s) at one state -- what the solver actually integrates"
function _sys1(T, ur, pphi, peta, PiB, tau, r, par = _PAR1)
    A = [zeros(5,5)]; S = zeros(5)
    Fluidum.matrix1d_visc!(A, S, [T, ur, pphi, peta, PiB], tau, (r,), par)
    A[1], S
end

@testset "1+1D viscous -- field layout and parities" begin
    f = Fluidum.viscous_1d()
    @test length(f.field) == 5
    @test collect(Fluidum.names(f)) == [:temperature, :ur, :piphiphi, :pietaeta, :piB]
    # r = 0 is a SYMMETRY AXIS, not a boundary: the radial velocity is odd there and every scalar
    # or diagonal-tensor channel is even. A parity error here is a kink in the first cells, which
    # is easy to mistake for a physical feature -- hence an explicit assertion.
    @test Fluidum.get_index(:pietaeta, f) == 4
    @test Fluidum.get_index(:ur, f)       == 2
    # the charm layout is this one plus (alpha, nu^r), in that order
    g = Fluidum.HQ_viscous_1d()
    @test collect(Fluidum.names(g))[1:5] == collect(Fluidum.names(f))
    @test collect(Fluidum.names(g))[6:7] == [:α, :nur]
end

@testset "1+1D viscous -- shear closure" begin
    # In radial symmetry pi^{mu nu} is diagonal in the parallel-transported triad (l, phihat,
    # etahat), so the two stored traceless channels fix the third: pi^l_l = -(pi^phi_phi+pi^eta_eta).
    # What is asserted is that the tensor those numbers describe is traceless and u-orthogonal.
    g = Diagonal([-1.0, 1.0, 1.0, 1.0])
    for ur in (-2.0, -0.5, 0.0, 0.7, 3.0), (pphi, peta) in ((0.1, -0.2), (-0.03, 0.07), (0.5, 0.5))
        ut  = sqrt(1 + ur^2); pll = -(pphi + peta)
        P   = [pll*ur^2  pll*ur*ut 0.0  0.0
               pll*ur*ut pll*ut^2  0.0  0.0
               0.0       0.0       pphi 0.0
               0.0       0.0       0.0  peta]
        u   = [ut, ur, 0.0, 0.0]
        sc  = max(maximum(abs, P), 1e-30)
        @test abs(tr(g*P))/sc                  < 1e-13
        @test maximum(abs, (u'*g*P)')/sc       < 1e-13
    end
end

@testset "1+1D viscous -- characteristic speeds against closed forms" begin
    # At rest, pi = Pi = 0: the largest Israel-Stewart characteristic is the MIS front
    #     v = sqrt(c_s^2 + (4/3) eta/(w tau_pi) + zeta/(w tau_Pi)),  w = e + p = T s.
    for T in (0.12, 0.2, 0.3, 0.45, 0.6)
        eta, taupi, zeta, tauPi = _coeffs1(T)
        Ar, _ = _sys1(T, 0.0, 0.0, 0.0, 0.0, 1.0, 5.0)
        lam   = maximum(abs, real.(eigvals(Ar)))
        v     = sqrt(_cs2(T, _EOS1) + (4/3)*eta/(_enth(T, _EOS1)*taupi) +
                     zeta/(_enth(T, _EOS1)*tauPi))
        @test isapprox(lam, v; rtol = 1e-12)
        @test v < 1.0                                   # and the front is subluminal
    end
    # With every dissipative coefficient off, the extreme characteristics are the boosted sound cone.
    for T in (0.15, 0.3, 0.5), ur in (-1.5, -0.4, 0.0, 0.9, 2.5)
        Ar, _ = _sys1(T, ur, 0.0, 0.0, 0.0, 1.0, 5.0, _PAR1I)
        lam   = sort(real.(eigvals(Ar)))
        v     = ur/sqrt(1 + ur^2); c = sqrt(_cs2(T, _EOS1))
        @test isapprox(lam[1],   (v - c)/(1 - v*c); atol = 1e-10)
        @test isapprox(lam[end], (v + c)/(1 + v*c); atol = 1e-10)
    end
end

@testset "1+1D viscous -- the derived matrix reproduces the shipped kernel" begin
    # `Julia/tools/derive_2p1d_viscous.wls` re-derives Israel-Stewart from scratch. Its 1-D routine
    # must reproduce the SHIPPED production kernel entry by entry -- that control is what makes the
    # 2-D half of the same derivation trustworthy where it DISAGREES with the shipped 2-D matrix.
    # both delta_pipi values (2026-10-13): dpp = 0 (FLUIDUM_DELTA_PIPI=0) and the default 4/3 tau_pi
    worst = 0.0
    for T in (0.12, 0.35, 0.55), ur in (-1.3, 0.0, 2.2),
        (pphi, peta, PiB) in ((0.0,0.0,0.0), (0.02,-0.05,0.01)), tau in (0.4, 7.0), r in (0.05, 12.0)
        eta, taupi, zeta, tauPi = _coeffs1(T)
        phi = [T, ur, pphi, peta, PiB]
        for dpp in (0.0, (4/3)*taupi)
            a = Fluidum.one_d_viscous_matrix(        phi, tau, r, _prs(T,_EOS1), _ent(T,_EOS1), _dst(T,_EOS1), zeta, eta, taupi, tauPi, dpp)
            b = Fluidum.one_d_viscous_matrix_derived(phi, tau, r, _prs(T,_EOS1), _ent(T,_EOS1), _dst(T,_EOS1), zeta, eta, taupi, tauPi, dpp)
            for (x, y) in zip(a, b)
                worst = max(worst, maximum(abs, Array(x) .- Array(y))/max(maximum(abs, x), 1e-30))
            end
        end
    end
    @test worst < 1e-12
end

@testset "1+1D viscous -- Bjorken is the 0+1D Israel-Stewart (+ delta_pipi) system" begin
    # The one configuration with a closed-form answer: no radial dependence, u^r = 0. The mixed
    # convention is pinned here -- pi^eta_eta -> -4 eta/(3 tau), pi^phi_phi -> +2 eta/(3 tau),
    # Pi -> -zeta/tau. On Bjorken theta = 1/tau, so the shear row the matrix encodes is
    #     tau_pi d_tau pi^eta_eta + pi^eta_eta + 4 eta/(3 tau) + delta_pipi pi^eta_eta/tau = 0
    # with delta_pipi = (4/3) tau_pi by default and 0 under FLUIDUM_DELTA_PIPI=0 (2026-10-13; until
    # then this kernel was plain IS). No tau_pipi, no lambda_piPi. Both switch values are asserted,
    # each against its own ODE; the two answers differ by ~20 % in pi^eta_eta on the state reached
    # below, so each is a sharp statement, not a smoke test.
    T0, tau0, tauf = 0.45, 0.4, 4.0
    function solve_ref(f)
        function rhs!(du, U, p, tau)
            T, peta, PiB = U
            eta, taupi, zeta, tauPi = _coeffs1(T)
            du[1] = (-(_ene(T,_EOS1) + _prs(T,_EOS1) + PiB + peta)/tau)/(T*_dst(T,_EOS1))
            du[2] = (-peta - 4eta/(3tau) - f*(4/3)*taupi*peta/tau)/taupi
            du[3] = (-PiB  - zeta/tau)/tauPi
        end
        solve(ODEProblem(rhs!, [T0, 0.0, 0.0], (tau0, tauf)), Tsit5(); reltol = 1e-12, abstol = 1e-14)
    end
    had = Fluidum.DELTA_PIPI[]
    try
        for fac in (1.0, 0.0)
            Fluidum.DELTA_PIPI[] = fac
            d   = DiscreteFields(Fluidum.viscous_1d(),
                                 CartesianDiscretization(Fluidum.OriginInterval(16, 20.0)), Float64)
            phi = set_array(r -> T0, :temperature, d)
            for fl in (:ur, :piphiphi, :pietaeta, :piB); set_array!(phi, r -> 0.0, fl, d); end
            sol = Fluidum.oneshoot(d, Fluidum.matrix1d_visc!, _PAR1, phi, (tau0, tauf);
                                   reltol = 1e-11, abstol = 1e-13)

            fin = sol[end]; k = size(fin, 2) ÷ 2
            Tr, petar, PiBr = solve_ref(fac)(tauf)
            @test isapprox(fin[1, k], Tr;    rtol = 1e-9)
            @test isapprox(fin[4, k], petar; rtol = 1e-8)
            @test isapprox(fin[5, k], PiBr;  rtol = 1e-8)
            # pi^phi_phi = -pi^eta_eta/2 must EMERGE from two independently evolved rows, not be imposed
            @test abs(fin[3, k] + fin[4, k]/2)/abs(fin[4, k]) < 1e-9
            # uniformity and rest are exact fixed points of the discretisation
            @test maximum(abs, fin[1, :] .- fin[1, k])/fin[1, k] < 1e-11
            @test maximum(abs, fin[2, :])                        < 1e-11
            # and the OTHER switch value's ODE is far away: the term is really there (or really absent)
            @test abs(fin[4, k]/solve_ref(1.0 - fac)(tauf)[2] - 1) > 0.05
        end
    finally
        Fluidum.DELTA_PIPI[] = had
    end
end

@testset "1+1D viscous -- the 1-D and 2-D transport paths agree" begin
    # `matrix1d_visc!` reaches the transport coefficients through a scalar entropy, `matrix2d_visc!`
    # through a whole `Thermodynamic`. They are the same formulas and must give the same numbers --
    # the two entry points once accepted DISJOINT viscosity models for exactly this reason
    # (FLUIDUM_2P1D_HQ.md section 5).
    for T in (0.12, 0.25, 0.4, 0.6)
        dpt  = Fluidum.pressure_derivative(T, Val(1), _EOS1)
        dptt = Fluidum.pressure_derivative(T, Val(2), _EOS1)
        th   = thermodynamic(T, _EOS1)
        @test isapprox(Fluidum.viscosity(T, dpt, _PAR1.shear),
                       Fluidum.viscosity(T, th,  _PAR1.shear); rtol = 1e-14)
        @test isapprox(Fluidum.τ_shear(T, dpt, _PAR1.shear),
                       Fluidum.τ_shear(T, th,  _PAR1.shear);  rtol = 1e-14)
        @test isapprox(Fluidum.bulk_viscosity(T, dpt, _PAR1.bulk),
                       Fluidum.bulk_viscosity(T, th,  _PAR1.bulk); rtol = 1e-14)
        @test isapprox(Fluidum.τ_bulk(T, dpt, dptt, _PAR1.bulk),
                       Fluidum.τ_bulk(T, th, _PAR1.bulk);          rtol = 1e-14)
    end
end
