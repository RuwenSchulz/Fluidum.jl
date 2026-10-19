# The conformal set and the two-way medium + charm system (src/Matrix/two_way.jl).
# Derivation: deriv/derive_two_way.wls.
using Fluidum, Test, StaticArrays, LinearAlgebra, OrdinaryDiffEq

@testset "two-way system" begin
    F = Fluidum
    CH = BoltzmannCharm()
    CEOS = ConformalEOS(47.5)
    PARC = F.FluidProperties(CEOS, QGPViscosity(0.16, 0.2), ZeroBulkViscosity(), HQdiffusion(0.116, 1.5))

    @testset "ConformalEOS" begin
        th = thermodynamic(0.3, CEOS); a = conformal_aSB(CEOS)
        @test th.pressure ≈ a * 0.3^4 && th.pressure_derivative[1] ≈ 4a * 0.3^3 && th.pressure_hessian[1] ≈ 12a * 0.3^2
        t2 = thermodynamic(0.3, -2.0, CEOS)
        @test t2.pressure == th.pressure && t2.pressure_derivative == (th.pressure_derivative[1], 0.0)
        @test t2.pressure_hessian == (th.pressure_hessian[1], 0.0, 0.0)
        @test F.pressure_derivative(0.3, Val(3), CEOS) ≈ 24a * 0.3
        @test thermodynamic(1, ConformalEOS(40)).pressure isa Float64
        @test 3 * F.pressure(0.3, CEOS) ≈ F.energy_density(0.3, -1.0, CEOS)       # ε = 3P
    end

    @testset "charm EOS" begin
        l = charm_quark_list()
        for (T, α) in ((0.3, -2.0), (0.16, 0.5), (0.5, -10.0))
            q1 = two_way_charm_thermo(T, α, CH, 1.5); q2 = two_way_charm_thermo(T, α, l, 1.5)
            @test all(isapprox.(values(q1), values(q2); rtol = 1e-13))
            tr = F.hq_transport(T, α, F.Heavy_Quark(l.particle_list, 1.0), HQdiffusion(0.116, 1.5))
            @test q1.h ≈ tr.h rtol = 1e-13                                     # τ_n = D_s h/T, same h
            @test F.tw_Ds(T, HQdiffusion(0.116, 1.5)) * q1.cQ / (2T) ≈ tr.tauM rtol = 1e-13
            δ = 1e-5; hp = (two_way_charm_thermo(T + δ, α, CH, 1.5).h - two_way_charm_thermo(T - δ, α, CH, 1.5).h) / 2δ
            @test q1.dhdT ≈ hp rtol = 1e-7
            @test F.tw_alpha_from_density(T, q1.n, CH) ≈ α atol = 1e-10
        end
    end

    @testset "conformal ≡ generic with ConformalEOS (1+1D S5)" begin
        for c in (false, true)
            ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 0.01 * 0.0, 1e-4, -5e-5, 2e-5)
            ϕg = SVector(ϕ[1], ϕ[2], ϕ[3], ϕ[4], 0.0, ϕ[5], ϕ[6], ϕ[7], ϕ[8], ϕ[9])
            rc = tw_system_1p1d(:two_way, ϕ, 1.3, 2.7, PARC, TwoWayOptions(; charm = CH, cM = c))
            rg = tw_system_1p1d(:two_way, ϕg, 1.3, 2.7, PARC, TwoWayOptions(; charm = CH, cM = c, generic = true))
            k = [1, 2, 3, 4, 6, 7, 8, 9, 10]
            for j in 1:2
                @test maximum(abs.(rg[j][k, k] .- rc[j])) <= 1e-13 * maximum(abs.(rc[j]))
            end
            @test maximum(abs.(rg[3][k] .- rc[3])) <= 1e-13 * maximum(abs.(rc[3]))
        end
    end

    @testset "c = 0: second moment passive; c = 1: live" begin
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 1e-4, 1e-4, -5e-5, 2e-5)
        r0 = tw_system_1p1d(:two_way, ϕ, 1.3, 2.7, PARC, TwoWayOptions(; charm = CH, cM = false))
        r1 = tw_system_1p1d(:two_way, ϕ, 1.3, 2.7, PARC, TwoWayOptions(; charm = CH, cM = true))
        @test all(iszero, r0[1][1:6, 7:9]) && all(iszero, r0[2][1:6, 7:9])
        @test any(!iszero, r1[1][1:6, 7:9]) || any(!iszero, r1[2][1:6, 7:9])
    end

    @testset "closures: zero allocations" begin
        m1 = two_way_1p1d!(; charm = CH)
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 0.01, 0.001, -0.0005, 0.0002)
        A = (zeros(MMatrix{9,9,Float64}),); S = zeros(MVector{9,Float64})
        m1(A, S, ϕ, 1.0, (3.0,), PARC); m1(A, S, ϕ, 1.0, (3.0,), PARC)
        @test (@allocated m1(A, S, ϕ, 1.0, (3.0,), PARC)) == 0
        @test all(isfinite, A[1]) && all(isfinite, S)
        m2 = two_way_2p1d!(; charm = CH)
        ϕ2 = SVector(0.3, 0.2, -0.1, 0.001, -0.002, 0.0003, -2.0, 0.01, -0.005, 0.001, -0.0005, 0.0002, 0.0001)
        A2 = (zeros(MMatrix{13,13,Float64}), zeros(MMatrix{13,13,Float64})); S2 = zeros(MVector{13,Float64})
        m2(A2, S2, ϕ2, 1.0, (3.0, 1.0), PARC); m2(A2, S2, ϕ2, 1.0, (3.0, 1.0), PARC)
        @test (@allocated m2(A2, S2, ϕ2, 1.0, (3.0, 1.0), PARC)) == 0
    end

    @testset "Bjorken" begin
        sol = two_way_bjorken(:ideal, [0.45], (0.4, 5.0); params = PARC)
        @test sol(5.0)[1] * 5.0^(1/3) ≈ 0.45 * 0.4^(1/3) rtol = 1e-9
        α0 = F.tw_alpha_from_density(0.45, 0.3, CH)
        for c in (false, true)
            s = two_way_bjorken(:two_way, [0.45, 0.0, α0, 0.0, 0.0], (0.4, 5.0); params = PARC, charm = CH, cM = c, dmn_eps = 0.0)
            u = s(5.0)
            @test 5.0 * two_way_charm_thermo(u[1], u[3], CH, 1.5).n ≈ 0.4 * 0.3 rtol = 1e-8
        end
    end

    @testset "cold edge: decoupling factor and two_way_problem" begin
        o = TwoWayOptions(; charm = CH)
        @test tw_edge_factor(0.2, o) == 1.0 && tw_edge_factor(0.156, o) == 1.0
        @test tw_edge_factor(0.12, o) == o.edge_floor && 0 < tw_edge_factor(0.146, o) < 1
        @test tw_edge_factor(0.12, TwoWayOptions(; charm = CH, edge = false)) == 1.0
        ϕ = SVector(0.3, 0.2, 0.001, -0.002, -2.0, 1e-4, 1e-4, -5e-5, 2e-5)      # above T_edge: untouched
        r0 = tw_system_1p1d(:two_way, ϕ, 1.3, 2.7, PARC, TwoWayOptions(; charm = CH, edge = false))
        r1 = tw_system_1p1d(:two_way, ϕ, 1.3, 2.7, PARC, o)
        @test all(r0[k] == r1[k] for k in 1:3)
        disc = CartesianDiscretization(F.OriginInterval(60, 12.0))
        d = DiscreteFields(two_way_conformal_1d(), disc, Float64)
        phi = set_array(r -> 0.12 + 0.25 / (1 + exp((r - 5) / 0.6)), :temperature, d)
        set_array!(phi, r -> F.tw_alpha_from_density(0.12 + 0.25 / (1 + exp((r - 5) / 0.6)), 0.2 * exp(-(r / 3)^2) + 1e-12, CH), :α, d)
        p1 = two_way_problem(d, PARC, phi, (0.4, 3.0); charm = CH)
        p0 = two_way_problem(d, PARC, phi, (0.4, 3.0); charm = CH, corona_flux = false)
        du1 = similar(phi); du0 = similar(phi); p1.f(du1, phi, nothing, 0.4); p0.f(du0, phi, nothing, 0.4)
        hot = [phi[1, i] >= o.T_edge for i in axes(phi, 2)]
        @test du1[:, hot] == du0[:, hot] && du1[[1, 2, 3, 4, 6, 7, 8, 9], :] == du0[[1, 2, 3, 4, 6, 7, 8, 9], :]
        @test du1[5, .!hot] != du0[5, .!hot]
        g = p1.f.f
        allocs(f, a...) = (f(a...); @allocated f(a...))        # in a function: no global-scope boxing
        @test allocs(F._tw_corona_1d!, du1, phi, 0.4, g) == 0
        sol = solve(p1, Tsit5(); reltol = 1e-8, abstol = 1e-10, dtmax = 0.01)
        @test string(sol.retcode) == "Success" && all(isfinite, sol.u[end])
    end

    @testset "field layouts" begin
        @test length(F.names(two_way_1d())) == 10 && length(F.names(two_way_conformal_1d())) == 9
        @test length(F.names(two_way_2d())) == 14 && length(F.names(two_way_conformal_2d())) == 13
        @test length(F.names(HQ_2p1d_BG_m2_7f_init())) == 7
    end
end
