# ==============================================================================================
# test/test_2p1d_viscous.jl -- the 2+1D viscous solver, `matrix2d_visc!`.
#
# Until 2026-09-03 the whole coverage of this solver was one `@test all(isfinite.(result))` in the
# "2d viscous" testset below, and that is how ten wrong entries in the shear rows survived. What is
# asserted here is what a smoke test cannot see: closed-form characteristic speeds, the isotropy the
# defect broke, and the Bjorken limit -- all of it cheap enough to run on every `Pkg.test()`.
#
# The heavy, quantitative version of the same statements (sound attenuation against the exact MIS
# dispersion relation, Gubser, 2-D against the 1-D cylindrical solver) is
#     Julia/Fluidum.jl/bench/bench_2p1d_viscous.jl
# and the derivation side is
#     Julia/tools/derive_2p1d_viscous.wls  +  Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl
# ==============================================================================================
using Test, LinearAlgebra
# no Printf: under `Pkg.test` the environment is Fluidum's deps plus [targets] test = ["Test"],
# and Printf is in neither, so importing it here turns the whole suite red while still looking
# green under `--project=Julia`, whose outer environment does have it.
using Fluidum, OrdinaryDiffEq
const _F = Fluidum

_prs(T, eos)  = Fluidum.pressure(T, eos)
_ent(T, eos)  = Fluidum.pressure_derivative(T, Val(1), eos)
_dst(T, eos)  = Fluidum.pressure_derivative(T, Val(2), eos)
_ene(T, eos)  = T*_ent(T, eos) - _prs(T, eos)
_enth(T, eos) = T*_ent(T, eos)
_cs2(T, eos)  = _ent(T, eos)/(T*_dst(T, eos))

const _EOS  = FluiduMEoS()
const _PAR  = Fluidum.FluidProperties(_EOS, QGPViscosity(0.1, 0.2),
                                      SimpleBulkViscosity(0.083, 15.0), ZeroDiffusion())
const _PARI = Fluidum.FluidProperties(_EOS, ZeroViscosity(), ZeroBulkViscosity(), ZeroDiffusion())

function _coeffs(T, par = _PAR)
    th = thermodynamic(T, par.eos)
    (Fluidum.viscosity(T, th, par.shear),      Fluidum.τ_shear(T, th, par.shear),
     Fluidum.bulk_viscosity(T, th, par.bulk),  Fluidum.τ_bulk(T, th, par.bulk))
end
function _sys(T, ux, uy, pyy, pzz, pxy, PiB, tau, par = _PAR)
    A = [zeros(7,7), zeros(7,7)]; S = zeros(7)
    Fluidum.matrix2d_visc!(A, S, [T,ux,uy,pyy,pzz,pxy,PiB], tau, (1.0,1.0), par)
    A[1], A[2], S
end
_spec(A) = sort(real.(eigvals(A)))

@testset "2+1D viscous -- field layout" begin
    f = Fluidum.viscous_2d()
    @test length(f.field) == 7
    @test collect(Fluidum.names(f)) == [:temperature, :ux, :uy, :piyy, :pizz, :pixy, :piB]
    @test Fluidum.get_index(:pizz, f) == 5      # the mixed pi^eta_eta -- see the Bjorken testset
end

@testset "2+1D viscous -- shear closure" begin
    # the three stored components must close onto a traceless, u-orthogonal pi^{mu nu} through
    #   pi^{xx} = [ -(pi^{yy}+pi^{zz})(u^t)^2 + 2 u^x u^y pi^{xy} + (u^y)^2 pi^{yy} ] / (1 + (u^y)^2)
    g = Diagonal([-1.0, 1.0, 1.0, 1.0])
    for ux in (-2.0, 0.0, 0.7), uy in (-1.5, 0.0, 2.0),
        (pyy, pzz, pxy) in ((0.1,-0.2,0.05), (-0.03,0.07,-0.11))
        ut  = sqrt(1 + ux^2 + uy^2)
        pxx = (-(pyy + pzz)*ut^2 + 2ux*uy*pxy + uy^2*pyy)/(1 + uy^2)
        ptx = (ux*pxx + uy*pxy)/ut; pty = (ux*pxy + uy*pyy)/ut
        ptt = (ux*ptx + uy*pty)/ut
        P = [ptt ptx pty 0.0; ptx pxx pxy 0.0; pty pxy pyy 0.0; 0.0 0.0 0.0 pzz]
        u = [ut, ux, uy, 0.0]; sc = maximum(abs, P)
        @test abs(tr(g*P))/sc < 1e-13
        @test maximum(abs, (u'*g*P)')/sc < 1e-13
    end
end

@testset "2+1D viscous -- characteristic speeds against closed forms" begin
    for T in (0.15, 0.30, 0.60)
        eta, taupi, zeta, tauPi = _coeffs(T); w = _enth(T, _EOS)
        # the Israel-Stewart front, at rest
        vf = sqrt(_cs2(T,_EOS) + (4/3)*eta/(w*taupi) + zeta/(w*tauPi))
        @test maximum(abs, _spec(_sys(T,0,0,0,0,0,0,2.0)[1])) ≈ vf rtol=1e-12
        # The transverse shear channel.  Israel-Stewart gives sqrt(eta/(w tau_pi)) = sqrt(C_s),
        # T- and EoS-independent.  This assertion used to demand sqrt(2 C_s) -- the closed form of
        # the WRONG equation: the legacy matrix's one-sided sigma^{xy} makes that channel a factor
        # sqrt(2) too fast along x and absent along y.  Corrected when the derived matrix became
        # the default on 2026-09-03; the legacy value is asserted below so the defect stays covered.
        @test minimum(abs.(_spec(_sys(T,0,0,0,0,0,0,2.0)[1]) .- sqrt(0.2))) < 1e-12
        let had = Fluidum.VISC_2D_DERIVED[]
            try
                Fluidum.VISC_2D_DERIVED[] = false
                @test minimum(abs.(_spec(_sys(T,0,0,0,0,0,0,2.0)[1]) .- sqrt(0.4))) < 1e-12
            finally
                Fluidum.VISC_2D_DERIVED[] = had
            end
        end
        # with viscosity off, the extreme characteristics are the relativistic sound cone
        for U in (0.0, 1.2)
            c = sqrt(_cs2(T,_EOS)); v = U/sqrt(1+U^2)
            ev = _spec(_sys(T,U,0,0,0,0,0,2.0,_PARI)[1])
            @test ev[end] ≈ (v+c)/(1+v*c) rtol=1e-11
            @test ev[1]   ≈ (v-c)/(1-v*c) rtol=1e-11
        end
        # A_t must be invertible everywhere the solver goes
        @test isfinite(sum(_sys(T,0.3,-0.2,0.01,-0.02,0.003,0.005,2.0)[1]))
    end
end

# ── THE TWO ASSERTIONS THAT WOULD HAVE CAUGHT THE pi^{xy} DEFECT FROM COLD ──────────────────────
# Both are sub-millisecond and need no PDE solve. They are stated against the DERIVED matrix, which
# is the statement about the physics; the shipped matrix fails both, which is why
# `FLUIDUM_2D_DERIVED=1` is recommended for new 2+1D work (see the note atop src/Matrix/2d_viscous.jl).
@testset "2+1D viscous -- isotropy and the symmetric sigma^{xy} (derived matrix)" begin
    had = Fluidum.VISC_2D_DERIVED[]
    try
        Fluidum.VISC_2D_DERIVED[] = true
        for T in (0.20, 0.45)
            Ax, Ay, _ = _sys(T, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
            # at rest and isotropic, the x- and y-flux spectra must coincide
            @test maximum(abs, _spec(Ax) .- _spec(Ay)) < 1e-11
            # and pi^{xy} must be driven by the SYMMETRIC gradient: eta/tau_pi on each half
            eta, taupi, _, _ = _coeffs(T)
            @test Ax[6,3] ≈ eta/taupi rtol=1e-10
            @test Ay[6,2] ≈ eta/taupi rtol=1e-10
        end
    finally
        Fluidum.VISC_2D_DERIVED[] = had
    end
end

@testset "2+1D viscous -- Bjorken is the 0+1D Israel-Stewart system" begin
    # A transversely uniform state has no transverse gradient, so the 2-D solver must integrate
    #   De = -(e + P + Pi + pi^eta_eta)/tau,  tau_pi D pi^eta_eta = -pi^eta_eta - 4 eta/(3 tau),
    #   tau_Pi D Pi = -Pi - zeta/tau
    # written here from the physics rather than from the matrix. This is what pins the conventions:
    # pi^{zz} is the MIXED pi^eta_eta = tau^2 pi^{eta eta}, and the relaxation is PLAIN
    # Israel-Stewart -- a delta_pipi term of FiVoHydro's size would show at the 10 % level.
    tau0, tau1, T0 = 0.4, 3.0, 0.35
    function rhs!(du, U, p, tau)
        T, pe, PB = U
        eta, taupi, zeta, tauPi = _coeffs(T)
        du[1] = (-(_ene(T,_EOS) + _prs(T,_EOS) + PB + pe)/tau)/(T*_dst(T,_EOS))
        du[2] = (-pe - 4eta/(3tau))/taupi
        du[3] = (-PB - zeta/tau)/tauPi
    end
    ref = solve(ODEProblem(rhs!, [T0,0.0,0.0], (tau0,tau1)), Tsit5(); reltol=1e-12, abstol=1e-14)(tau1)

    d = DiscreteFields(Fluidum.viscous_2d(),
        CartesianDiscretization(Fluidum.SymmetricInterval(24, 12.0),
                                Fluidum.SymmetricInterval(24, 12.0)), Float64)
    phi = set_array((x,y)->T0, :temperature, d)
    for f in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    u = Fluidum.oneshoot(d, Fluidum.matrix2d_visc!, _PAR, phi, (tau0,tau1);
                         reltol=1e-11, abstol=1e-13)(tau1)
    ic = 13
    @test u[1,ic,ic] ≈ ref[1] rtol=1e-10
    @test u[5,ic,ic] ≈ ref[2] rtol=1e-10
    @test u[7,ic,ic] ≈ ref[3] rtol=1e-10
    # tracelessness is DYNAMICAL: pi^{yy} and pi^{zz} are advanced by separate rows
    @test u[4,ic,ic] ≈ -u[5,ic,ic]/2 rtol=1e-11
    # a uniform state must stay uniform and at rest
    @test maximum(abs, u[1,2:25,2:25] .- u[1,ic,ic])/u[1,ic,ic] < 1e-13
    @test maximum(abs, u[2,2:25,2:25]) < 1e-13
    @test maximum(abs, u[3,2:25,2:25]) < 1e-13

    # and with viscosity off, tau*s is conserved on the Bjorken solution
    phi = set_array((x,y)->T0, :temperature, d)
    for f in (:ux,:uy,:piyy,:pizz,:pixy,:piB); set_array!(phi, (x,y)->0.0, f, d); end
    ui = Fluidum.oneshoot(d, Fluidum.matrix2d_visc!, _PARI, phi, (tau0,tau1);
                          reltol=1e-11, abstol=1e-13)(tau1)
    @test tau1*_ent(ui[1,ic,ic], _EOS) ≈ tau0*_ent(T0, _EOS) rtol=1e-9
end

@testset "2+1D viscous -- the 1-D and 2-D transport paths agree" begin
    # These were DISJOINT method sets until 2026-09-02: `matrix2d_visc!` passes a whole
    # `Thermodynamic` and `matrix1d_visc!` a scalar entropy, and each entry point worked with
    # exactly the viscosity model the other could not take. The formulas are identical; this is
    # what makes sure they stay identical.
    for T in (0.15, 0.25, 0.35, 0.5)
        th = thermodynamic(T, _EOS); dpt = _ent(T,_EOS); dptt = _dst(T,_EOS)
        @test Fluidum.viscosity(T, dpt, _PAR.shear)     ≈ Fluidum.viscosity(T, th, _PAR.shear)     rtol=1e-12
        @test Fluidum.bulk_viscosity(T, dpt, _PAR.bulk) ≈ Fluidum.bulk_viscosity(T, th, _PAR.bulk) rtol=1e-12
    end
end
