# Fluidum

[![Build Status](https://github.com/fafafrens/Fluidum.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fafafrens/Fluidum.jl/actions/workflows/CI.yml?query=branch%3Amain)

To install the package and its dependencies one can use 
```julia
import Pkg
Pkg.add(url = "https://github.com/fafafrens/Fluidum.jl")
```
After this one can use the package simply by using 
```julia
using Fluidum
```

## Example 1+1 dimension
A simple example of running heavy quark simulations 
(actually is not simple at all since it is not understandable: what Observable does?
what runFluidum ? 
Why we need 2 equations of state ? 
)

```julia
using Fluidum
#define equation of state and transport coefficients
eos = Heavy_Quark()
viscosity=QGPViscosity(0.2,0.2)
bulk=SimpleBulkViscosity(0.1,15.0)
diffusion=ZeroDiffusion()
#collect eos and transport coefficient in the FluidProperties struct
params=Fluidum.FluidProperties(eos,viscosity,bulk,diffusion)
#define a radial grid
gridpoints=100
rmax=30
disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
#define the name and number of fields, with appropriate parity
oned_visc_hydro = oned_visc_hydro=Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB)
    )
    
disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
#define a example function for the temperature
function temperature(r)
       0.4*1/(exp(r/7)+1 )+0.0001
end
#set the temperature field
phi=set_array((x)->temperature(x),:temperature,disc_fields); 
#define a time range for the evolution
tspan = (0.2,1000)
#evolve until freeze-out temperature Tfo
freeze_out_routine(disc_fields,Fluidum.matrix1d_visc!,params,phi,tspan;Tfo=0.1565)
```

## Example 2+1 dimension 
Instead you can also solve in 2 dimension 
```julia
using Fluidum
#set up the fluid poroperties 
fluidpropery=FluidProperties(FluiduMEoS(),SimpleShearViscosity(0.1,.1),ZeroBulkViscosity(),ZeroDiffusion())

# the convention here are T, ux, uy, piyy, pizz, pixy, piB this has to match with the matrix defined
twod_visc_hydro=Fields(
NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
NDField((:ghost,:ghost),(:ghost,:ghost),:piB)
)

#we define a 2 cartesian grid form -25 to 25 50 point each dimension 
discretization=CartesianDiscretization(Fluidum.SymmetricInterval(50,25.),Fluidum.SymmetricInterval(50,25.))

# we prepare the field with the discretization
twod_visc_hydro_discrete=DiscreteFields(twod_visc_hydro,discretization,Float64)

 #we define some random intial condition 
function temperature(r)
       0.4(1+0.3rand())/(exp(abs(r)-10)+1 )+0.01
end
#we set the array corresponding to the temperature 
phi=set_array((x,y)->temperature(hypot(x,y)),:temperature,twod_visc_hydro_discrete);

using Plots
#here to undertand whath we are doing 
plot(phi[1,:,25])

#this is the time span 
tspan=(0.4,20)


# this create the ODEProblem  and feed it diffenential equations 
res=oneshoot(twod_visc_hydro_discrete,Fluidum.matrix2d_visc!,fluidpropery,phi,tspan)

#plot the solution 
using Plots
plot(phi[1,:,25])
plot!(res[end][1,:,25])
```

## The 1+1D and 2+1D viscous and charm solvers

Boost-invariant Milne solvers, cylindrical (1+1D) and transverse Cartesian (2+1D). Documentation:

| | |
|---|---|
| `examples/1p1d_viscous/README.md` | how to run the 1+1D viscous and charm solvers, the mixed `pi^eta_eta` convention, the axis parity, and three worked examples |
| `examples/2p1d_viscous/README.md` | the same for the 2+1D transverse solver |
| `bench/README.md` | the two closed-form benchmarks: 12 checks in 1+1D, 19 in 2+1D |
| `deriv/` | the Wolfram scripts that generate the `*_generated.jl` matrices |

The drivers:

| driver | fields | what it is |
|---|---|---|
| `matrix1d_visc!` | `(T, u^r, pi^phi_phi, pi^eta_eta, Pi)` | the medium in cylindrical Milne coordinates |
| `matrix1d_visc_HQ!` | the five above + `(alpha, nu^r)` | 1+1D medium and charm, charm passive |
| `matrix2d_visc!` | `(T, u^x, u^y, pi^{yy}, pi^eta_eta, pi^{xy}, Pi)` | the medium, Israel-Stewart |
| `matrix2d_HQ_BG!` | `(alpha, nu^x, nu^y)` | the charm first moment on a prescribed background |
| `matrix2d_visc_HQ_BG!` | the seven above + `(alpha, nu^x, nu^y)` | medium and charm together, charm passive |
| `matrix2d_HQ_BG_m2!` | `(alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy}, piQ^eta_eta, Pi_Q)` | charm first and second moment on a prescribed background |

Two defaults to know before comparing with other codes:

* `FLUIDUM_2D_DERIVED` is on by default: the 2+1D viscous kernels use the matrix generated by
  `deriv/derive_2p1d_viscous.wls`. The legacy hand-written matrix (`FLUIDUM_2D_DERIVED=0`) has ten
  wrong entries in the two shear rows; use it only to reproduce older results.
* `matrix2d_HQ_BG_m2!` takes `consistent_fm` and defaults it to `true` (the full-grad-P charm
  first moment). The other two 2+1D charm drivers use the grad-alpha row and have no switch. State
  the choice explicitly in any cross-code run.

## Citing 
When using Fluidum.jl for research, teaching or similar, please cite our work, see [CITATIONS.bib](CITATIONS.bib).
