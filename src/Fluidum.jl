module Fluidum

using ForwardDiff:Dual,value,partials
using Bessels
using ForwardDiff
using DifferentiationInterface
using DelimitedFiles: readdlm
using DelimitedFiles
using StructArrays
using IntervalSets
using MuladdMacro
using StaticArrays
using Artifacts 

#pkg used for intial conditions


#pkg used by the fluid evolution. 
using Interpolations
using OrdinaryDiffEq
using Statistics
using JLD2


using LinearAlgebra
using LoopVectorization
using Surrogates
using StatsBase

using FastChebInterp


using TensorCast
using SimpleNonlinearSolve

using QuadGK
# NumericalIntegration is declared in Project.toml and referenced by `multiplicity_analytic`
# (fluidevo/spectra.jl) as `NumericalIntegration.integrate`, but was never imported — that call
# was an undefined GlobalRef. (2026-09-02)
using NumericalIntegration

using HCubature
using UnPack



const root_particle_lists=artifact"particle_lists"


const root_kernels=artifact"kernels"

# ── Process-wide charm-sector switches ───────────────────────────────────────────────────────
# Every one of them is a Ref populated in __init__, NOT `const X = parse(..., ENV[...])`. A plain
# const is evaluated at PRECOMPILE time and baked into the .ji, so later processes silently reuse
# whatever value was current when the image was built and the env var has no effect whatsoever —
# verified the hard way: two solves at HQ_TAUN_SCALE 1 and 6 came out bit-identical
# (max|Δν^r| = 0.0). __init__ runs on every module load, so the Refs track the environment of the
# process. Read ONCE per process (2026-10-13; until then four of them were `get(ENV, …)` inside the
# per-cell RHS): set the env var before `using Fluidum`, or assign the Ref in-process
# (`Fluidum.HQ_CM_BACKREACTION[] = true`) before building the matrix. Setting ENV after the module
# is loaded does nothing.
#
# HQ_TAUN_SCALE: diagnostic multiplier on the heavy-quark diffusion relaxation time τ_n, applied by
# `hq_taun` / `hq_transport` (EquationOfState/heavyquark.jl) and therefore by EVERY charm matrix
# that builds τ_n. Default 1.0 = BARE = production (since 2026-07-21 τ_diffusion_hadron itself is
# bare, Eq. (30) of 2205.07692). Set HQ_TAUN_SCALE=1/6 ONLY to reproduce the historical ÷g_hq value
# for the convention study (diag_taun_hydro_vs_langevin.jl). No production path sets it.
const HQ_TAUN_SCALE = Ref(1.0)
# FLUIDUM_HQ_MOMENTUM_DIMS: 3 (production) or 2 = charm τ_n/τ_M matched in the 2-D momentum measure of
# a transverse-only Langevin ensemble — see `hq_transport` in EquationOfState/heavyquark.jl. Diagnostic.
const HQ_MOMENTUM_DIMS = Ref(3)
# HQ_CM_BACKREACTION: the 1st↔2nd-moment back-coupling c_M. false (DEFAULT, env "0" or unset) ⇒
# c_M = 0, which is what every existing product uses; true (env "1") ⇒ c_M = D_s/T, the elliptic
# branch (FLUIDUM_ATTRACTOR_AUDIT.md §3). Until 2026-10-13 the default was c_M = D_s/T and only the
# literal "0" gave zero, so every production driver had to pin it.
const HQ_CM_BACKREACTION = Ref(false)
# HQC_CM_NFLOOR: the density floor that damps the c_M force in the charm-free tail
# (HQ_const_BG_consistent.jl). Irrelevant at c_M = 0. Default 1e-6 fm^-3.
const HQC_CM_NFLOOR = Ref(1e-6)
# LP1_M2_SOURCES: the source content of the consistent second moment's rows 3-5
# (HQ_const_BG_consistent_m2.jl): :consistent (default), :legacy or :nabla_nu.
const HQ_M2_SOURCES = Ref(:consistent)
# FIVO_IS2_TAUM_BARE: the diagnostic `_rel` τ_M/η_M branch of the RETRACTED HQ_const_BG_2nd_moment.jl.
const HQ_TAUM_BARE_DIAG = Ref(false)
# FLUIDUM_ALLOW_RETRACTED: lets a RETRACTED matrix run (to reproduce an old artifact). Default false.
const ALLOW_RETRACTED = Ref(false)
# FLUIDUM_DELTA_PIPI: the DNMR delta_pipi of the 5-field 1-D medium (`one_d_viscous_matrix`), its
# generated twin and the 2-D medium (`two_d_viscous_matrix_derived`), as a factor on (4/3) tau_pi
# (`delta_pipi`, Matrix/1d_viscous_cilindrical.jl). 1.0 (DEFAULT, env unset or "1") = 4/3 tau_pi,
# RS 2026-10-13 ("everywhere"; the 7-field :MIS kernel always had it); 0.0 (env "0") = plain IS,
# which reproduces every pre-2026-10-13 result of those kernels bit for bit.
const DELTA_PIPI = Ref(1.0)
# FLUIDUM_HQ_M2_MIXTURE: the consistent charm second moment's EOS coefficients (τ_M ∝ C, d ln C/d ln T,
# A = I41/P0, d ln I31/d ln T) for a MULTI-SPECIES charm list. true (DEFAULT, env unset or "1") = the
# list's own species sums (`hq_m2_mixture`, EquationOfState/heavyquark.jl), shared by the m2 drivers
# AND the two-way wrappers; false (env "0") = the pre-2026-10-13 single-quark Bessel forms at
# m = diffusion.mass in the drivers and C = h K4K2/K3², A = m² + 5Th in the wrappers, bit for bit.
# A ONE-SPECIES list (unit charge, mass = diffusion.mass — LP1/O+O's charm.data) and BoltzmannCharm
# take the old code path either way, so they are bit-identical at both values.
const HQ_M2_MIXTURE = Ref(true)

function _env_bool(name, default::String)
    s = get(ENV, name, default)
    s == "0" && return false
    s == "1" && return true
    error("$name must be \"0\" or \"1\", got \"$s\"")
end

function __init__()
    HQ_TAUN_SCALE[] = parse(Float64, get(ENV, "HQ_TAUN_SCALE", "1.0"))
    HQ_MOMENTUM_DIMS[] = parse(Int, get(ENV, "FLUIDUM_HQ_MOMENTUM_DIMS", "3"))
    HQ_MOMENTUM_DIMS[] in (2, 3) || error("FLUIDUM_HQ_MOMENTUM_DIMS must be 2 or 3")
    HQ_MOMENTUM_DIMS[] == 3 || @warn "Fluidum: charm τ_n/τ_M matched in the 2-D momentum measure — DIAGNOSTIC MODE, not production"
    HQ_TAUN_SCALE[] == 1.0 || @warn "Fluidum: τ_n scaled by HQ_TAUN_SCALE — DIAGNOSTIC MODE, not production" scale=HQ_TAUN_SCALE[]
    HQ_CM_BACKREACTION[] = _env_bool("HQ_CM_BACKREACTION", "0")
    HQ_CM_BACKREACTION[] && @warn "Fluidum: c_M = D_s/T (HQ_CM_BACKREACTION=1) — the elliptic branch; every existing product is c_M = 0"
    HQC_CM_NFLOOR[] = parse(Float64, get(ENV, "HQC_CM_NFLOOR", "1e-6"))
    HQ_M2_SOURCES[] = _parse_m2_sources(get(ENV, "LP1_M2_SOURCES", "consistent"))
    HQ_TAUM_BARE_DIAG[] = _env_bool("FIVO_IS2_TAUM_BARE", "0")
    ALLOW_RETRACTED[] = _env_bool("FLUIDUM_ALLOW_RETRACTED", "0")
    ALLOW_RETRACTED[] && @warn "Fluidum: FLUIDUM_ALLOW_RETRACTED=1 — retracted matrices will run (reproduction of old artifacts only)"
    # FLUIDUM_2D_DERIVED: the re-derived 2+1D viscous matrix, DEFAULT ON since 2026-09-03 -- ten
    # entries of the legacy one are wrong (Matrix/2d_viscous.jl).  Set 0 to get the legacy matrix,
    # which is what everything produced before that date used.  Same reason as above for reading it
    # HERE and not at top level.
    VISC_2D_DERIVED[] = get(ENV, "FLUIDUM_2D_DERIVED", "1") != "0"
    VISC_2D_DERIVED[] || @warn "Fluidum: 2+1D viscous hydro on the LEGACY matrix (FLUIDUM_2D_DERIVED=0) — ten entries are wrong; this is for reproducing pre-2026-09-03 results only"
    DELTA_PIPI[] = _env_bool("FLUIDUM_DELTA_PIPI", "1") ? 1.0 : 0.0
    DELTA_PIPI[] == 1.0 || @warn "Fluidum: delta_pipi = 0 in the 5-field 1-D and the 2-D viscous medium (FLUIDUM_DELTA_PIPI=0) — plain Israel-Stewart, for reproducing pre-2026-10-13 results only; the 7-field :MIS kernel keeps 4/3 tau_pi"
    HQ_M2_MIXTURE[] = _env_bool("FLUIDUM_HQ_M2_MIXTURE", "1")
    HQ_M2_MIXTURE[] || @warn "Fluidum: charm second moment with SINGLE-QUARK coefficients for multi-species lists (FLUIDUM_HQ_M2_MIXTURE=0) — for reproducing pre-2026-10-13 list results only"
end

const fmGeV= 1/0.1973261
const invfmGeV= 1/fmGeV
const invfmGeV3=(1/fmGeV)^3
const fmGeV3=fmGeV^3


include("EquationOfState/thermodynamic.jl")
include("EquationOfState/thermodynamic_perturbation.jl")
include("EquationOfState/EquationofStatetype.jl")
include("EquationOfState/TransportCoefficienttype.jl")


include("EquationOfState/IdealQCD.jl")
include("EquationOfState/FluiduMEoS.jl")
include("EquationOfState/conformalEOS.jl")

include("EquationOfState/charm_HRG.jl")
include("EquationOfState/simpletransport.jl")
include("EquationOfState/Analytic.jl")

include("EquationOfState/heavyquark.jl")


include("fluidevo/detector_struct.jl")
include("fluidevo/inverse_function.jl")



include("fluidevo/discretization.jl") 
include("fluidevo/picewisefunction.jl") 
include("fluidevo/discreteFields.jl")  

include("fluidevo/generating_kernels_3d.jl")
#include("fluidevo/generating_kernels.jl")
include("fluidevo/particle_dictionary.jl")


include("fluidevo/spectra_fastreso_dict_HQ.jl")
include("fluidevo/spectra.jl")

include("fluidevo/test_functions.jl")
include("fluidevo/map_profile.jl")
include("fluidevo/initial_fields.jl")

include("Matrix/1d_viscous_HQ_cilindrical_fugacity.jl")
include("Matrix/1d_viscous_HQ_cilindrical_gamma.jl")

include("Matrix/1d_ideal_HQ_cilindrical_density.jl")
include("Matrix/1d_ideal_HQ_cilindrical_fugacity.jl")
include("Matrix/1d_ideal_HQ_cilindrical_p_br_density.jl")
include("Matrix/1d_ideal_HQ_cilindrical_p_br_fugacity.jl")
include("Matrix/HQ_const_BG.jl")
# ⛔ RETRACTED 2026-10-13 (rows 3-5 wrong) — kept for reproduction; refuses to run unless forced.
include("Matrix/HQ_const_BG_2nd_moment.jl")
# The consistent 1-D charm sector, in the package build since 2026-10-13 (until then drivers
# `Base.include`d them at run time, which still works: the files define only functions). Order
# matters for reading only: _consistent_m2 calls into _consistent, which calls the shipped
# second-moment rows of the file above for its :legacy control.
include("Matrix/HQ_const_BG_consistent.jl")
include("Matrix/HQ_const_BG_consistent_m2.jl")
# 2+1D charm on a prescribed background.  The matrices are GENERATED — regenerate with
#   wolframscript -file Julia/tools/derive_hq_2p1d.wls
include("Matrix/HQ_2p1d_BG_generated.jl")
# the CONSISTENT (full-grad-P) 2-D first moment -- derive_hq_consistent_2p1d.wls
include("Matrix/HQ_2p1d_BG_consistent_generated.jl")
include("Matrix/HQ_2p1d_BG.jl")
# the 2-D consistent second moment (needs hqc_enthalpy_per_particle from HQ_const_BG_consistent.jl
# and hq2d_matrices[_consistent] from the two generated files above)
include("Matrix/HQ_2p1d_BG_m2.jl")

include("Matrix/1d_viscous_cilindrical.jl")
# 2026-09-03: the from-scratch Israel-Stewart derivation (Julia/tools/derive_2p1d_viscous.wls).
# It reproduces `one_d_viscous_matrix` entry by entry, which is what licenses it to be used as the
# reference for the 2-D matrix in Projects/FluidumValidation/bench_fluidum_2p1d_viscous.jl.
include("Matrix/viscous_generated.jl")

include("Matrix/2d_viscous.jl")
include("Matrix/2d_viscous_fugacity.jl")

# 2026-10-13 (P3b): the conformal set and the TWO-WAY medium + charm system in 0+1D, 1+1D, 2+1D.
# The matrices are GENERATED — regenerate with  wolframscript -file Julia/tools/derive_two_way.wls
include("Matrix/two_way_generated.jl")
include("Matrix/two_way.jl")

include("wraps2.jl")
#include("plotting_routines.jl")

const  detector_collection=(;ALICE=detector("ALICE",6.62,7.00,0.0757,"Pb_Pb"),
RHIC =detector("RHIC" ,7.,4.23,0.005968,"Au_Au"),
ALICE1 =detector("ALICE1",     6.62 ,    7.00 ,        0.0463,	"Pb_Pb")
)

const detector_dict=Dict(
:ALICE_NNPDF=>detector(:ALICE,6.62,7.00,0.0757,:Pb_Pb),
:ALICE_CTEQ =>detector(:ALICE,6.62,7.00 ,0.0463,:Pb_Pb),
:ALICE_average =>detector(:ALICE,6.62,7.00, 0.061,:Pb_Pb),
:RHIC =>detector(:RHIC ,7.,4.23,0.005968,:Au_Au),
)

 
export detector_collection,detector_dict

export NDField, Fields, OriginInterval, CartesianDiscretization, DiscreteFields
# 2026-09-02: seven names were exported that are defined NOWHERE in the package —
#   spectra_analitic (typo for spectra_analytic), initialize_fields_free_HQ, map_initial_profile,
#   Profiles2, free_charm (renamed to free_hadron in 471c28c, Apr 2026), InverseFuction,
#   local_qnm_spectrum_BG_second_moment (deleted with 218 lines of HQ_const_BG_2nd_moment.jl in
#   030fcc9, Mar 2026).
# `using Fluidum` bound all seven to nothing, and `free_charm` is what aborted test/runtests.jl at
# the FIRST testset — testsets 2-7 ("fluid_properties", "Cheb", "initialconditions", "2d viscous",
# "causality check", "artifacts") had therefore never run on this branch.  Removed; the real
# `spectra_analytic` / `multiplicity_analytic` are exported in their place.
export set_array, set_array!, freeze_out_routine, fo_integral, jgemvavx!, oneshoot, oneshoot_debug, test_integral_cauchy, SplineInterp, spectra_analytic, multiplicity_analytic
export spectra, spectra_lf, multiplicity, multiplicity_lf, DiscreteFields, TabulatedData, initialize_fields, FreezeOutResult, FreezeOutResultPerturbation
export get_profile, Profiles #RunFluidum_hf, RunFluidum_array, save_to_h5, RunFluidum_lf, SetFluidProperties, FluidParameters

"""

    thermodynamic(T[,μ],x::EquationofState)
Give back the pressure, the pressure gradient and pressure hessian at given T and if specified at (T,μ) for given equation of state
"""
function thermodynamic
end

"""
    pressure(T[,μ],x::Equation of state)
Give back the pressure at given T and if specified at (T,μ) for given equation of state


"""
function pressure
end

"""
    pressure_derivative(T[,μ],Val(N)[,Val(M)],x::EquationofState)
Return the pressure gradient/hessian

# Arguments
- `T::Number`: Temperature
- `μ::Number`: Chemical potential
- `Val(N::Integer)`: N is the order of derivatives with respect to T ranging from 0 to 2
- `Val(M::Integer)`: M is the order of derivatives with respect to μ ranging from 0 to 2
- `x::EquationOfState`: equation of state

# Examples
Compute ∂_μ p for LatticeQCD at (T=0.5,μ=0.3)
```jldoctest
>julia LQCD=LatticeQCD()
LatticeQCD()
>julia pressure_derivative(0.5,0.3,Val(0),Val(1),LQCD)
0.024387044999999035
```
"""
function pressure_derivative
end

"""
    OneDPicewiseEquationOfState(eos1::EquationofState, int1::AbstractInterval,eos2::EquationOfState,int2::AbstractInterval)
Give back p_1/2(T), if T is in interval int1/2, else give back 0

"""
function OneDPicewiseEquationOfState
end

"""
    TwoDPicewiseEquationOfState(eos1::EquationofState, int1::AbstractInterval,eos2::EquationOfState,int2::AbstractInterval)
Give back p_1/2(T,μ), if (T,μ) is in interval int1/2, else give back 0
"""
function TwoDPicewiseEquationOfState
end

"""
    viscosity(T[,μ][,x::EquationOfState,th::Thermodynamic],y::TransportProperty)
Give back shear viscosity η/s and relaxation time τ_s 

# Arguments
- `T::Number`: Temperature
- `μ::Number`: Chemical potential
- `x::EquationOfState`: equation of state. Either this or thermodynamic has to be specified.
- `th::Thermodynamic`: thermodynamic. Either this or the equation of state has to be specified.
- `y::TransportProperty`: transport TransportProperty

# Examples
Compute the shear viscosity for T=0.3 GeV and μ=0.2 GeV for a hadron resonance gas
```jldoctest
>julia hrg=HadronResonaceGas()
HadronResonaceGas()
>julia vis=EquationsOfStates.SimpleShearViscosity(0.16,0.2)
SimpleShearViscosity{Float64}(0.16, 0.2)
>julia th=thermodynamic(0.3,0.2,hrg)
p(T,μ)=0.06585468437527325
∇p(T,μ)=(1.581024911648665, 0.07161301647333912)
∇²p(T,μ)=(30.990803012108447, 1.3978520817058184, 30.990803012108447)
>julia viscosity(0.3,0.2,hrg,vis)
0.0499163967709561
>julia viscosity(0.3,0.2,th,vis)
0.0499163967709561
```
"""
function viscosity
end

"""
    IdealQCD(Nf,Nc)
Set up fluid properties using the pressure of an ideal QCD gas with Nf flavours and Nc colors
"""
function IdealQCD
end

"""
    FluiduMEoS
Set up fluid properties using the pressure from a fit to HRG at low temperature and lattice QCD at high temperature
"""
function FluiduMEoS
end

"""
    HadronResonaceGas([,name_file,Maxmass,Minmass,condition])
Set up fluid properties using the pressure of a hadron resonance gas

# Arguments
- `name_file::String`: Name of file with particle information. Standart value is "EquationofState/particles.data"
- `Maxmass::Number`: Mass cut-off in GeV, heavier particles will be excluded. Standart value is 2.1 
- `Minmass::Number`: Mass cut-off in GeV, lighter particles will be excluded. Standart value is 0.05 
- `condition`: Custom condition to exclude particles
"""
function HadronResonaceGas
end


"""
    waleckacondition
Condition to exclude particles from hadron resonance gas already included in walecka model (ρ,ω,p^+,p^-)
"""
function waleckacondition
end

"""
    LatticeQCD
Set up fluid properties using the pressure from lattice QCD (paper ref)
"""
function LatticeQCD
end

"""
    WaleckaModel
Set up the fluid properties using the pressure from the waleckamodel
"""
function WaleckaModel
end

"""
    HadronResonaceGasNew
Set up the fluid properties from a hadron resonance gas based on the therminator particle file
"""
function HadronResonaceGasNew
end

"""
    SimpleBulkViscosity(ζs,Cζ)
Set up the bulk viscosity and corresponding relaxation time

# Arguments
- `ζs::Number`: Dimensionless value of maximum of bulk Bulk viscosity (ζ/s)_max, recommended value is 0.05
- `Cζ::Number`: Dimensionless value for magnitude of τ_bulk (∝ 1/C), recommended value is 15
"""
function SimpleBulkViscosity
end


"""
    SimpleShearViscosity(ηs,Cη)
Set up the shear viscosity and corresponding relaxation time

# Arguments
- `ηs::Number`: Value of shear viscosity over entropy ration η/s


TBW
"""
function SimpleShearViscosity
end

export thermodynamic,pressure, pressure_derivative , TwoDPicewiseEquationOfState ,OneDPicewiseEquationOfState, energy_density
export FluidProperties, viscosity, τ_shear,bulk_viscosity,τ_bulk,diffusion,τ_diffusion
export IdealQCD, FluiduMEoS, HadronResonaceGas,waleckacondition,LatticeQCD, HadronResonaceGasNew, HadronResonaceGas_ccbar
export SimpleBulkViscosity,SimpleShearViscosity,SimpleDiffusionCoefficient , ZeroViscosity, ZeroBulkViscosity
export Analytic, Gluing,Thermodynamic, EquationOfState, readresonancelist
export Heavy_Quark, HQdiffusion, free_hadron, QGPViscosity,ZeroDiffusion
export matrix2d_HQ_BG!, HQ_2p1d_BG_init, hq2d_matrices, axisymmetric_bg_2d
export matrix2d_visc_HQ_BG!, HQ_viscous_2d, hq2d_live_rows
# The 2-D second moment's field layout and driver, and the consistent 1-D closures. The drivers
# were not exported until 2026-10-13 because their files were not part of the build.
export HQ_2p1d_BG_m2_init, matrix2d_HQ_BG_m2!
export matrix1d_visc_HQ_BG_consistent!, matrix1d_visc_HQ_BG_consistent_5f!, matrix1d_visc_HQ_BG_consistent_m2_5f!
# 2026-10-13 (P3b): the conformal set + the two-way system (Matrix/two_way.jl), their layouts and EOS
export ConformalEOS, conformal_aSB, BoltzmannCharm, charm_quark_list, two_way_charm_thermo, TwoWayOptions
export two_way_1p1d!, tw_ideal_1p1d!, tw_viscous_1p1d!, tw_hq_fm_bg_1p1d!, tw_hq_m2_bg_1p1d!
export two_way_2p1d!, tw_ideal_2p1d!, tw_viscous_2p1d!, tw_hq_fm_bg_2p1d!, tw_hq_m2_bg_2p1d!
export tw_system_0p1d, tw_system_1p1d, tw_system_2p1d, tw_bjorken_rhs, two_way_bjorken, tw_alpha_from_density
export two_way_problem, tw_edge_factor
export tw_ideal_1d, tw_viscous_conformal_1d, two_way_1d, two_way_conformal_1d
export tw_ideal_2d, tw_viscous_conformal_2d, HQ_2p1d_BG_m2_7f_init, two_way_2d, two_way_conformal_2d
# type piracy

Bessels.besselk0(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(Bessels.besselk0(value(d)), -Bessels.besselk1(value(d)) * partials(d))

Bessels.besselk1(d::Dual{T,V,N}) where {T,V,N} = Dual{T}(Bessels.besselk1(value(d)),-( Bessels.besselk0(value(d))+Bessels.besselk1(value(d))/value(d)) * partials(d))

end
