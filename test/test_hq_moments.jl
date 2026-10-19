# The charm second moment and the first-moment rows it uses: the retracted-matrix switch, the c_M
# default, the shared transport helper, passivity, the 2-D -> 1-D reduction, zero allocations.
using Fluidum, Test, StaticArrays, LinearAlgebra, QuadGK

const HQM_CHARM_FILE = let p = joinpath(mktempdir(), "charm.data")
    write(p, "# Name Mass Gamma deg Spin Isospin I3 Nq Ns Naq Nas Nc Nac MC\n" *
             "Charm 1.5 0.0 6 0.5 0. 0. 0. 0. 0. 0. 1. 0. 4\n")
    p
end
const HQM_EOS = Fluidum.Heavy_Quark(readresonancelist(name_file = HQM_CHARM_FILE), 1.0)
const HQM_PARAMS = Fluidum.FluidProperties(HQM_EOS, QGPViscosity(0.1, 0.2),
                                           SimpleBulkViscosity(0.083, 15.0), HQdiffusion(0.11634, 1.5))
const HQM_BG1 = (tau, r) -> (0.31, 0.42, -0.055, -0.021, 0.13, 0.07)
const HQM_BG2 = (tau, x, y) -> (0.3, 0.4, -0.2, -0.05, -0.02, 0.01, 0.05, 0.1, 0.02, -0.03, 0.04, 0.08)
const HQM_PHI5 = SVector(-2.0, 0.031, 1.0e-3, -4.0e-4, 2.0e-4)

"Run a 1-D closure at (tau, r); return (A_t^-1 A_x, A_t^-1 src)."
function hqm_eval1(m, ϕ; tau = 2.3, r = 3.1, params = HQM_PARAMS)
    N = length(ϕ)
    A = (zeros(MMatrix{N,N,Float64}),); S = zeros(MVector{N,Float64})
    m(A, S, ϕ, tau, (r,), params)
    return SMatrix(A[1]), SVector(S)
end
function hqm_eval2(m, ϕ; tau = 1.7, x = 3.0, y = -1.0, params = HQM_PARAMS)
    A = (zeros(MMatrix{8,8,Float64}), zeros(MMatrix{8,8,Float64})); S = zeros(MVector{8,Float64})
    m(A, S, ϕ, tau, (x, y), params)
    return SMatrix(A[1]), SMatrix(A[2]), SVector(S)
end
function hqm_alloc(m, ϕ, X, N, ndim)
    A = ntuple(_ -> zeros(MMatrix{N,N,Float64}), ndim); S = zeros(MVector{N,Float64})
    m(A, S, ϕ, 1.2, X, HQM_PARAMS); m(A, S, ϕ, 1.2, X, HQM_PARAMS)
    return @allocated m(A, S, ϕ, 1.2, X, HQM_PARAMS)
end
"Run `f` with a Ref temporarily set."
function hqm_with(f, ref::Ref, v)
    old = ref[]; ref[] = v
    try
        return f()
    finally
        ref[] = old
    end
end

@testset "HQ second moment" begin

@testset "retracted shipped second moment refuses to run" begin
    hqm_with(Fluidum.ALLOW_RETRACTED, false) do
        @test_throws ErrorException Fluidum.matrix1d_visc_HQ_BG_second_moment!(; background_fields = HQM_BG1)
        A = (zeros(MMatrix{5,5,Float64}),); S = zeros(MVector{5,Float64})
        @test_throws ErrorException Fluidum.matrix1d_visc_HQ_BG_second_moment!(
            A, S, HQM_PHI5, 2.3, (3.1,), HQM_PARAMS; background_fields = HQM_BG1)
        m = Fluidum.matrix1d_visc_HQ_BG_second_moment!(; background_fields = HQM_BG1, allow_retracted = true)
        Ai, Si = hqm_eval1(m, HQM_PHI5)
        @test all(isfinite, Ai) && all(isfinite, Si)
    end
    hqm_with(Fluidum.ALLOW_RETRACTED, true) do      # same as FLUIDUM_ALLOW_RETRACTED=1
        @test Fluidum.matrix1d_visc_HQ_BG_second_moment!(; background_fields = HQM_BG1) isa Function
    end
end

@testset "c_M defaults to 0" begin
    if get(ENV, "HQ_CM_BACKREACTION", "0") == "0"
        @test Fluidum.HQ_CM_BACKREACTION[] == false
    end
    @test Fluidum._env_bool("HQ_CM_BACKREACTION_UNSET_IN_ANY_ENV", "0") == false
    withenv("HQM_TEST_FLAG" => "1") do
        @test Fluidum._env_bool("HQM_TEST_FLAG", "0") == true
    end
    withenv("HQM_TEST_FLAG" => "yes") do
        @test_throws ErrorException Fluidum._env_bool("HQM_TEST_FLAG", "0")
    end
    # with c_M = 0 the row-2 c_M entries are exactly zero; with c_M = D_s/T they are not
    hqm_with(Fluidum.HQ_CM_BACKREACTION, false) do
        b = Fluidum.hqc_first_moment_block(HQM_PHI5, 2.3, 3.1, HQM_PARAMS, HQM_BG1(2.3, 3.1), 1e-6, 200.0, false)
        @test b.atc == 0 && b.axc == 0 && b.csrc == 0
    end
    hqm_with(Fluidum.HQ_CM_BACKREACTION, true) do
        b = Fluidum.hqc_first_moment_block(HQM_PHI5, 2.3, 3.1, HQM_PARAMS, HQM_BG1(2.3, 3.1), 1e-6, 200.0, false)
        @test b.atc != 0 && b.axc != 0 && b.csrc != 0
    end
end

@testset "passivity at c_M = 0: m2 fields do not feed the first moment" begin
    hqm_with(Fluidum.HQ_CM_BACKREACTION, false) do
        ϕb = SVector(HQM_PHI5[1], HQM_PHI5[2], -3.0e-3, 7.0e-4, 1.1e-3)    # same (α, ν), other moments
        for mk in (() -> Fluidum.matrix1d_visc_HQ_BG_consistent_5f!(; background_fields = HQM_BG1),
                   () -> Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1),
                   () -> Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1, sources = :legacy),
                   () -> Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1, sources = :nabla_nu))
            A, S = hqm_eval1(mk(), HQM_PHI5)
            Ab, Sb = hqm_eval1(mk(), ϕb)
            # after A_t^-1 (a dense 5x5 LU, which may pivot across the blocks) the decoupling
            # holds to roundoff, not bit for bit
            @test maximum(abs, A[1:2, 3:5]) <= 1e-14 * maximum(abs, A)
            @test isapprox(S[1:2], Sb[1:2]; rtol = 1e-13)
            @test isapprox(A[1:2, 1:2], Ab[1:2, 1:2]; rtol = 1e-13)
        end
        # on the assembled blocks, bit for bit: rows 1-2 carry no m2 column and no m2 field
        bg = HQM_BG1(2.3, 3.1)
        ϕb = SVector(HQM_PHI5[1], HQM_PHI5[2], -3.0e-3, 7.0e-4, 1.1e-3)
        b  = Fluidum.hqc_first_moment_block(HQM_PHI5, 2.3, 3.1, HQM_PARAMS, bg, 1e-6, 200.0, false)
        bb = Fluidum.hqc_first_moment_block(ϕb, 2.3, 3.1, HQM_PARAMS, bg, 1e-6, 200.0, false)
        @test b.atc == 0 && b.axc == 0 && b.csrc == 0
        @test b.At2 == bb.At2 && b.Ax2 == bb.Ax2 && b.src2 == bb.src2
        # 2-D: rows 1-3 never see columns 4-8
        ϕ8  = SVector(-2.0, 0.01, -0.02, 1e-3, 2e-4, -5e-4, -5e-4, 3e-4)
        ϕ8b = SVector(-2.0, 0.01, -0.02, -2e-3, 6e-4, 1e-3, 1e-3, -7e-4)
        for vt in (false, true)
            m = Fluidum.matrix2d_HQ_BG_m2!(; background_fields = HQM_BG2, vorticity = vt)
            Ax, Ay, S = hqm_eval2(m, ϕ8); _, _, Sb = hqm_eval2(m, ϕ8b)
            @test all(iszero, Ax[1:3, 4:8]) && all(iszero, Ay[1:3, 4:8])
            @test S[1:3] == Sb[1:3]
        end
    end
end

@testset "consistent_m2 :legacy reproduces consistent_5f bit for bit" begin
    for cm in (false, true), (tau, r, ϕ) in ((2.3, 3.1, HQM_PHI5), (0.6, 0.4, SVector(1.2, -0.04, 5e-3, 2e-3, -1e-3)),
                                              (4.7, 9.0, SVector(-4.5, 0.002, -1e-4, 3e-4, 6e-5)))
        hqm_with(Fluidum.HQ_CM_BACKREACTION, cm) do
            A5, S5 = hqm_eval1(Fluidum.matrix1d_visc_HQ_BG_consistent_5f!(; background_fields = HQM_BG1), ϕ; tau, r)
            Am, Sm = hqm_eval1(Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1,
                                                                            sources = :legacy), ϕ; tau, r)
            @test A5[1:2, :] == Am[1:2, :] && S5[1:2] == Sm[1:2]
            @test A5 == Am && S5 == Sm
        end
    end
    @test_throws ErrorException Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1, sources = :legacyy)
end

@testset "one (τ_n, τ_M, η_M, h) for every matrix" begin
    T, α, eos, dif = 0.27, -1.3, HQM_EOS, HQM_PARAMS.diffusion
    tr = Fluidum.hq_transport(T, α, eos, dif)
    z = dif.mass / T
    # the defaults reproduce the closed-form expressions bit for bit
    @test tr.taun == Fluidum.τ_diffusion_hadron(T, α, eos, dif) * 1.0
    @test tr.tauM == tr.taun * Fluidum.besselk(4, z) * Fluidum.besselk(2, z) / (2 * Fluidum.besselk(3, z)^2)
    @test tr.etaM == T * tr.taun / 2
    @test tr.h == tr.taun * T / tr.Ds
    @test Fluidum.hq_taun(T, α, eos, dif) == tr.taun
    @test Fluidum.hqc_enthalpy_per_particle(T, α, eos, dif) == (tr.h, tr.taun, tr.Ds)
    # both knobs act on every path identically
    for (ref, v) in ((Fluidum.HQ_MOMENTUM_DIMS, 2), (Fluidum.HQ_TAUN_SCALE, 0.5))
        hqm_with(ref, v) do
            t2 = Fluidum.hq_transport(T, α, eos, dif)
            @test t2.taun != tr.taun
            @test Fluidum.hq_taun(T, α, eos, dif) == t2.taun
            @test Fluidum.hqc_enthalpy_per_particle(T, α, eos, dif)[2] == t2.taun
            @test t2.h == t2.taun * T / t2.Ds              # tau_n = D_s h / T under every knob
            @test t2.etaM == T * t2.taun / 2
            # the 1-D first moment reads the same τ_n
            A0, _ = hqm_with(ref, ref === Fluidum.HQ_MOMENTUM_DIMS ? 3 : 1.0) do
                hqm_eval1(Fluidum.matrix1d_visc_HQ_BG!(; background_fields = HQM_BG1), SVector(α, 0.01))
            end
            A1, _ = hqm_eval1(Fluidum.matrix1d_visc_HQ_BG!(; background_fields = HQM_BG1), SVector(α, 0.01))
            @test A0 != A1
        end
    end
    # MOMENTUM_DIMS=2 rescales τ_M by (K3/K4)/λ2 in every matrix
    hqm_with(Fluidum.HQ_MOMENTUM_DIMS, 2) do
        t2 = Fluidum.hq_transport(T, α, eos, dif)
        λ1 = z * (z + 1) / (z^2 + 3z + 3); λ2 = z * (z^2 + 3z + 3) / (z^3 + 6z^2 + 15z + 15)
        K2, K3, K4 = Fluidum.besselk(2, z), Fluidum.besselk(3, z), Fluidum.besselk(4, z)
        @test t2.taun ≈ tr.taun * (K2 / K3) / λ1
        @test t2.tauM ≈ tr.tauM * (K3 / K4) / λ2
    end
end

@testset "2-D m2 rows reduce to the 1-D consistent m2 rows (axisymmetric)" begin
    # On the +x axis with u^y = 0 the 2-D m2 rates must reproduce the 1-D rates.
    STATES = [
        (tau=2.3, r=3.1, ur= 0.42, T=0.31, dtT=-0.055, drT=-0.021, drur= 0.13, dtur= 0.07,
         n=0.9, nu= 0.031, Ds=0.44, h=2.20, hp=-0.9,
         pl= 1.0e-3, pphi=-4.0e-4, PiQ= 2.0e-4, dta= 0.031, dra=-0.017, dtn= 2.3e-4, drn=-1.1e-4),
        (tau=1.1, r=0.7, ur=-0.30, T=0.45, dtT=-0.120, drT= 0.033, drur=-0.21, dtur= 0.15,
         n=1.4, nu=-0.070, Ds=0.30, h=1.90, hp=-1.4,
         pl=-7.0e-4, pphi= 3.0e-4, PiQ=-1.5e-4, dta=-0.020, dra= 0.040, dtn=-1.0e-4, drn= 3.0e-4),
        (tau=5.0, r=8.0, ur= 0.90, T=0.20, dtT=-0.010, drT=-0.005, drur= 0.02, dtur= 0.01,
         n=0.3, nu= 0.010, Ds=0.60, h=2.60, hp=-0.5,
         pl= 5.0e-4, pphi= 6.0e-4, PiQ= 9.0e-5, dta= 0.011, dra=-0.003, dtn= 5.0e-5, drn=-2.0e-5),
        (tau=0.6, r=1.2, ur=-1.30, T=0.55, dtT= 0.300, drT=-0.440, drur= 0.55, dtur=-0.62,
         n=2.2, nu= 0.190, Ds=0.12, h=1.75, hp=-2.2,
         pl= 9.0e-4, pphi= 1.0e-4, PiQ=-5.0e-4, dta= 0.090, dra=-0.060, dtn= 8.0e-4, drn=-4.0e-4),
    ]
    worst = 0.0
    for p in STATES, pj in (true, false)
        taun = p.Ds * p.h / p.T
        z = 1.5 / p.T
        tauM = taun * Fluidum.besselk(4, z) * Fluidum.besselk(2, z) / (2 * Fluidum.besselk(3, z)^2)
        etaM = p.T * taun / 2
        ux = p.ur; uT = sqrt(1 + ux^2)
        pxx = p.pl * uT^2                       # invert the l-hat projection on the +x axis
        pet = Fluidum.hq2d_piQ_eta(ux, 0.0, uT, pxx, 0.0, p.pphi)
        got = Fluidum.hq2d_m2_rows((0.0, p.nu, 0.0, pxx, 0.0, p.pphi, pet, p.PiQ), p.tau, p.T, ux, 0.0,
                                   p.dtT, p.drT, 0.0, p.dtur, p.drur, 0.0, 0.0, 0.0, p.ur / p.r,
                                   p.dra, 0.0, p.dta, p.dtn, p.drn, 0.0, 0.0, 0.0, p.nu / p.r,
                                   p.n, taun, p.Ds, p.h, p.hp, tauM, etaM, 1.5; projector = pj)
        at3, ax3, src3 = Fluidum.hqc_m2_rows((0.0, p.nu, p.pl, p.pphi, p.PiQ), p.tau, p.r, p.ur, p.T,
                                             p.dtT, p.drT, p.drur, p.dtur, p.n, taun, p.Ds, p.h, p.hp,
                                             tauM, etaM, 1.5; sources = :consistent)
        d1 = ntuple(i -> -(at3[i,1]*p.dta + at3[i,2]*p.dtn + ax3[i,1]*p.dra + ax3[i,2]*p.drn + src3[i]) /
                         at3[i, 2+i], 3)
        # the projection moves with the flow: dp_l/dτ = ∂_τπ^xx/u^τ² − 2π^xx u^r ∂_τu^r/u^τ⁴, and the
        # 2-D rate carries no advection, so the 1-D state's ∂_r p_l = 0 needs its −(u^x/u^τ)∂_xπ^xx back
        rate_xx = got[1] - (p.ur / uT) * (2 * pxx * p.ur * p.drur / uT^2)
        pl_2d = rate_xx / uT^2 - 2 * pxx * p.ur * p.dtur / uT^4
        rel(a, b) = abs(a - b) / max(abs(b), 1e-300)
        if pj      # the projector is part of the equation: the reduction holds only with it ON
            worst = max(worst, rel(pl_2d, d1[1]), rel(got[3], d1[2]), rel(got[4], d1[3]))
        else       # without it p_l must disagree
            @test rel(pl_2d, d1[1]) > 1e-8
        end
    end
    @test worst < 1e-10
end

@testset "zero allocations per cell" begin
    ϕ5 = HQM_PHI5; ϕ8 = SVector(-2.0, 0.01, -0.02, 1e-3, 2e-4, -5e-4, -5e-4, 3e-4)
    @test hqm_alloc(Fluidum.matrix1d_visc_HQ_BG!(; background_fields = HQM_BG1), SVector(-2.0, 0.01), (3.0,), 2, 1) == 0
    @test hqm_alloc(Fluidum.matrix1d_visc_HQ_BG_consistent!(; background_fields = HQM_BG1), SVector(-2.0, 0.01), (3.0,), 2, 1) == 0
    @test hqm_alloc(Fluidum.matrix1d_visc_HQ_BG_consistent_5f!(; background_fields = HQM_BG1), ϕ5, (3.0,), 5, 1) == 0
    for src in (:consistent, :legacy, :nabla_nu)
        @test hqm_alloc(Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1, sources = src),
                        ϕ5, (3.0,), 5, 1) == 0
    end
    for cf in (true, false), vt in (false, true)
        @test hqm_alloc(Fluidum.matrix2d_HQ_BG_m2!(; background_fields = HQM_BG2, consistent_fm = cf, vorticity = vt),
                        ϕ8, (3.0, -1.0), 8, 2) == 0
    end
end

@testset "multi-species list: second-moment coefficients = momentum quadrature, drivers ≡ wrappers" begin
    # three unit-charge species: the mixture values are unambiguous (no charge convention)
    L3 = Fluidum.HadronResonaceGas_ccbar(Fluidum.StructArray([
             Fluidum.Particle("D", 1.8646, 0.0, 4, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 421),
             Fluidum.Particle("Dst", 2.007, 0.0, 12, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 423),
             Fluidum.Particle("Lc", 2.2865, 0.0, 4, 0.5, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4122)]), 1.0)
    E3 = Fluidum.Heavy_Quark(mass = 1.5, hadron_list = L3)
    P3 = Fluidum.FluidProperties(E3, QGPViscosity(0.1, 0.2), SimpleBulkViscosity(0.083, 15.0), HQdiffusion(0.11634, 1.5))
    # I_nq by quadrature; ∂_T I_nq|_α = I_{n+1,q}/T² exactly
    function Inq(n, q, T, α, m, g)
        z = m / T; tm = acosh(1 + 900 / z)
        f(t) = m * cosh(t) * (m * sinh(t))^(2 + 2q) * (m * cosh(t))^(n - 2q - 1) * exp(-z * (cosh(t) - 1))
        return g * exp(α - z) / (2π^2 * prod(1:2:(2q+1))) * quadgk(f, 0.0, tm / 4, tm / 2, tm; rtol = 1e-14)[1]
    end
    worst = 0.0
    for T in (0.16, 0.25, 0.4), α in (-3.0, 1.0)
        S(n, q) = sum(Inq(n, q, T, α, p.Mass, p.Degeneracy) for p in L3.particle_list)
        h = S(3, 1) / S(2, 1); dh = (S(4, 1) - h * S(3, 1)) / T^2 / S(2, 1)
        C = S(5, 2) / S(4, 2); dC = (S(6, 2) - C * S(5, 2)) / T^2 / S(4, 2)
        ABr = T * S(4, 1) / (T^2 * S(3, 1))            # ∂ln I31/∂ln T, I31 = Σ I31_i
        A = S(4, 1) / S(2, 1)
        x = Fluidum.hq_m2_mixture(T, α, L3, 1.5)
        worst = max(worst, abs(x.h / h - 1), abs(x.dhdT / dh - 1), abs(x.C / C - 1), abs(x.dCdT / dC - 1),
                    abs(x.ABr / ABr - 1), abs((x.m2bar + 5T * x.h) / A - 1), abs(x.aI31 - 1))
        q = two_way_charm_thermo(T, α, E3, 1.5)
        worst = max(worst, abs(q.cQ / C - 1), abs(q.dcQdT / dC - 1), abs(q.aQ / A - 1), abs(T * q.dI31dT / q.I31 / ABr - 1))
        tr = Fluidum.hq_transport(T, α, E3, P3.diffusion)
        worst = max(worst, abs(tr.tauM * x.Cratio / (tr.Ds * C / (2T)) - 1))
    end
    @test worst < 1e-8
    # the drivers (their own FD h') ≡ the two-way S4 wrappers on that list, 1+1D and 2+1D
    hqm_with(Fluidum.HQ_CM_BACKREACTION, false) do
        X = SVector(-1.5, 0.004, 2e-3, -1e-3, 5e-4)
        A = (zeros(MMatrix{5,5,Float64}),); Sv = zeros(MVector{5,Float64})
        Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; dmn_eps = 0.0, background_fields = HQM_BG1)(A, Sv, X, 2.3, (3.1,), P3)
        B = (zeros(MMatrix{5,5,Float64}),); Q = zeros(MVector{5,Float64})
        tw_hq_m2_bg_1p1d!(; charm = E3, cM = false, dmn_eps = 0.0, nu_cap = 0.0, background_fields = HQM_BG1)(B, Q, X, 2.3, (3.1,), P3)
        @test maximum(abs.(B[1] .- A[1])) / maximum(abs.(A[1])) < 1e-6
        @test maximum(abs.(Q .- Sv)) / maximum(abs.(Sv)) < 1e-6
        # with the old coefficients the same comparison misses by percent
        Ao = (zeros(MMatrix{5,5,Float64}),); So = zeros(MVector{5,Float64})
        Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; dmn_eps = 0.0, background_fields = HQM_BG1, mixture = false)(Ao, So, X, 2.3, (3.1,), P3)
        @test maximum(abs.(So .- Sv)) / maximum(abs.(Sv)) > 1e-3
    end
    # one species (charm.data): bit-identical with the switch on and off
    for mx in (true, false)
        A = (zeros(MMatrix{5,5,Float64}),); Sv = zeros(MVector{5,Float64})
        Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1, mixture = mx)(A, Sv, HQM_PHI5, 2.3, (3.1,), HQM_PARAMS)
        A8 = (zeros(MMatrix{8,8,Float64}), zeros(MMatrix{8,8,Float64})); S8 = zeros(MVector{8,Float64})
        Fluidum.matrix2d_HQ_BG_m2!(; background_fields = HQM_BG2, mixture = mx)(A8, S8, SVector(-2.0, 0.01, -0.02, 1e-3, 2e-4, -5e-4, -5e-4, 3e-4), 1.7, (3.0, -1.0), HQM_PARAMS)
        if mx
            global HQM_ONE = (Array(A[1]), Array(Sv), Array(A8[1]), Array(A8[2]), Array(S8))
        else
            @test HQM_ONE == (Array(A[1]), Array(Sv), Array(A8[1]), Array(A8[2]), Array(S8))
        end
    end
    # zero allocations with a multi-species list
    @test (Fluidum.hq_m2_mixture(0.3, -1.0, L3, 1.5); @allocated Fluidum.hq_m2_mixture(0.3, -1.0, L3, 1.5)) == 0
    let m = Fluidum.matrix1d_visc_HQ_BG_consistent_m2_5f!(; background_fields = HQM_BG1)
        A = (zeros(MMatrix{5,5,Float64}),); Sv = zeros(MVector{5,Float64})
        m(A, Sv, HQM_PHI5, 1.2, (3.0,), P3); m(A, Sv, HQM_PHI5, 1.2, (3.0,), P3)
        @test (@allocated m(A, Sv, HQM_PHI5, 1.2, (3.0,), P3)) == 0
    end
    let m = Fluidum.matrix2d_HQ_BG_m2!(; background_fields = HQM_BG2, vorticity = true)
        ϕ8 = SVector(-2.0, 0.01, -0.02, 1e-3, 2e-4, -5e-4, -5e-4, 3e-4)
        A = (zeros(MMatrix{8,8,Float64}), zeros(MMatrix{8,8,Float64})); Sv = zeros(MVector{8,Float64})
        m(A, Sv, ϕ8, 1.2, (3.0, -1.0), P3); m(A, Sv, ϕ8, 1.2, (3.0, -1.0), P3)
        @test (@allocated m(A, Sv, ϕ8, 1.2, (3.0, -1.0), P3)) == 0
    end
end

end
