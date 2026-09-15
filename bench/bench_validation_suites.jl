#!/usr/bin/env julia
# ════════════════════════════════════════════════════════════════════════════════════════════════
# bench/bench_validation_suites.jl — RUN THE FOUR VALIDATION SUITES THAT NOTHING ELSE RUNS.
#
# WHY THIS FILE IS A SHIM AND NOT THE SUITES THEMSELVES
# -----------------------------------------------------
# Fluidum's deepest validation does not live in this package. It lives in
# `Julia/Projects/FluidumValidation/`, four files that between them cover every matrix entry point,
# the derivation of both viscous matrices, the charm sector and the 2+1D HQ system. And
# `run_suite.jl` discovers entries by walking `Julia/*/examples` and `Julia/*/bench` — so it cannot
# see them, and **nothing ran them but a person**. The repository's own rule is that a gate nobody
# runs is a comment.
#
# Moving them here would churn 23 references across seventeen files, two of which are historical
# logs this repository deliberately does not rewrite. So instead this file — which IS in a
# discovered `bench/` — invokes them as subprocesses and aggregates. Same pattern as
# `FiVoHydro.jl/test/run2d_gates.jl`, and for the same reason: the things being run cannot share a
# session (each sets its own module-level state and two of them load matrices INTO `Fluidum`).
#
# ⚠ THIS ONLY BECAME POSSIBLE ON 2026-09-14. Two of the four are AUDITS whose findings are failures:
# `bench_fluidum_full.jl` fails 7 gates and `bench_fluidum_attractor.jl` 4, all of them recorded
# findings from 2026-09-02 about entry points nothing uses. They therefore exited 1 unconditionally
# and could never have been wired into anything. They now DECLARE their failure sets
# (`KNOWN_FAILURES`) and exit 0 while the set is exactly the declared one — and 1 when it moves in
# either direction, a new failure OR a declared one starting to pass. That is what makes them
# runnable here.
#
#   julia --project=Julia Julia/Fluidum.jl/bench/bench_validation_suites.jl
#   FLUIDUM_AUDIT_FULL=1 …   also the heavy opt-in sections: the 2-D vs 1-D solve comparison and
#                            charm-number conservation (FLUIDUM_2D_FULL), and the MIS robustness
#                            map and freeze-out (FLUIDUM_BENCH_RUNS). Adds roughly 25 minutes.
# ════════════════════════════════════════════════════════════════════════════════════════════════
using Printf

const HERE = @__DIR__
const REPO = normpath(joinpath(HERE, "..", "..", ".."))
const VDIR = joinpath(REPO, "Julia", "Projects", "FluidumValidation")
const PROJ = joinpath(REPO, "Julia")
const FULL = get(ENV, "FLUIDUM_AUDIT_FULL", "0") == "1"

#        file                             what it covers                                     extra env when FULL
const SUITES = [
    ("bench_fluidum_full.jl",             "every matrix entry point + semi-analytic refs",    Dict("FLUIDUM_BENCH_RUNS" => "1")),
    ("bench_fluidum_attractor.jl",        "the charm sector on a closed-form flowing bg",     Dict{String,String}()),
    ("bench_fluidum_2p1d.jl",             "the 2+1D charm/HQ system",                         Dict("FLUIDUM_2D_FULL" => "1")),
    ("bench_fluidum_2p1d_viscous.jl",     "the 2+1D viscous DERIVATION and its audit",        Dict{String,String}()),
]

println("="^100)
println("Fluidum — the validation suites in Julia/Projects/FluidumValidation/   (",
        FULL ? "FULL" : "fast", " mode)")
println("`run_suite.jl` cannot discover them; this shim is how they get run.")
println("="^100)

function main()
    nfail = 0; nrun = 0
    for (file, what, extra) in SUITES
        path = joinpath(VDIR, file)
        if !isfile(path)
            # A listed suite whose file is gone is a FAILURE, not a skip — a silent skip is how a
            # validation suite stops existing without anyone noticing.
            @printf("  %-34s  %s\n", file, "FAIL (file missing)"); nfail += 1; continue
        end
        nrun += 1
        env = copy(ENV)
        FULL && for (k, v) in extra; env[k] = v; end
        io = IOBuffer(); t0 = time()
        ok = try
            run(pipeline(setenv(`$(Base.julia_cmd()) --project=$PROJ $path`, env);
                         stdout = io, stderr = io))
            true
        catch; false end
        out = String(take!(io))
        ok || (nfail += 1)
        @printf("  %-34s  %s  (%5.0f s)   %s\n", file, ok ? "PASS" : "FAIL", time() - t0, what)
        # echo the one line that says what the run concluded
        for ln in split(out, '\n')
            (occursin("passed,", ln) || occursin("gates passed", ln) ||
             occursin("failure set", ln) || occursin("KNOWN failures", ln)) &&
                println("        ", strip(ln))
        end
        ok || println(out)
    end
    println("="^100)
    @printf("%d/%d validation suites passed  (%s)\n", nrun - nfail, nrun, FULL ? "full" : "fast")
    exit(nfail == 0 ? 0 : 1)
end
main()
