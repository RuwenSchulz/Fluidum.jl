# Tests of src/attractor (this branch only). Run after runtests.jl:
#     julia --project test/attractor_tests.jl

using Fluidum, Test

@testset "attractor branch" begin
    @testset "free_hadron" begin
        @test isfinite(free_hadron(1, 0, 1, 1; m = 1.5)[1])
        @test all(isfinite, free_hadron(0.3, 0.0, 6, 1; m = 1.5))
    end

    @testset "HQdiffusion freezes the linear law at T_fo" begin
        @test HQdiffusion(0.116, 1.5) isa ConstDiffusion
        d = HQdiffusion(1.765, -0.159, 1.5)
        @test d.Tc == 0.156
        @test Fluidum.DsT(d, 0.10) ≈ 1.765 * 0.156 - 0.159
    end

    @testset "charm transport helpers" begin
        hq = Fluidum.Heavy_Quark(hadron_list = charm_quark_list())
        for dif in (ConstDiffusion(0.11634, 1.5), HQdiffusion(1.765, -0.159, 1.5)), T in (0.16, 0.3)
            tr = Fluidum.hq_transport(T, 0.2, hq, dif)
            @test tr.taun ≈ Fluidum.hq_taun(T, 0.2, hq, dif)
            @test tr.taun ≈ tr.Ds * tr.h / T
            @test tr.etaM ≈ T * tr.taun / 2
        end
    end

    @testset "conformal EoS with a fugacity slot" begin
        eos = ConformalEOS()
        th = Fluidum.thermodynamic(0.3, 0.7, eos)
        @test th.pressure ≈ Fluidum.thermodynamic(0.3, eos).pressure
        @test th.pressure_derivative[2] == 0
        @test all(iszero, th.pressure_hessian[2:3])
    end

    @testset "dict initializer converts the earlier second-moment keys" begin
        grid = (rmax = 10.0, gridpoints = 40)
        base = Dict{Symbol,Function}(:α => r -> 0.3 - 0.01r, :nu => r -> 1e-3r)
        old = merge(base, Dict{Symbol,Function}(:PiQr => r -> 1e-4, :PiQperp => r -> 2e-4, :PiQ => r -> -1e-4))
        new = merge(base, Dict{Symbol,Function}(:piphiphi_Q => r -> 2e-4, :pietaeta_Q => r -> -3e-4,
                                                :piB_Q => r -> -1e-4))
        lay = Fluidum.HQ_const_BG_second_moment_init
        @test Fluidum.initialize_fields(old, lay, grid).initial_field ≈
              Fluidum.initialize_fields(new, lay, grid).initial_field
        # `:nu` is the earlier key of the field `nur`
        renamed = merge(Dict{Symbol,Function}(:α => base[:α], :nur => base[:nu]),
                        Dict{Symbol,Function}(k => new[k] for k in (:piphiphi_Q, :pietaeta_Q, :piB_Q)))
        @test Fluidum.initialize_fields(new, lay, grid).initial_field == Fluidum.initialize_fields(renamed, lay, grid).initial_field
        @test any(!iszero, Fluidum.initialize_fields(new, lay, grid).initial_field)
    end
end
