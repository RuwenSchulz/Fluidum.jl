struct DiscreteInitialFields{A,B,C}    
    discretization::A
    discrete_field::B
    initial_field::C
end

"""
    HQ_viscous_1d()

Field layout for the 7-field 1+1D viscous hydro + charm system (`matrix1d_visc_HQ!`):
`(T, u^r, pi^phi_phi, pi^eta_eta, Pi, alpha, nu^r)` in cylindrical Milne. [`viscous_1d`](@ref)
plus the charm fugacity and its diffusion current; the Cartesian sibling is [`HQ_viscous_2d`](@ref).

`u^r` and `nu^r` are odd at the origin and everything else even, since r = 0 is a symmetry axis.
Charm is a passive tracer: it does not feed `T^{munu}`.
"""
@inline function HQ_viscous_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:α), 
    NDField((:odd,),(:ghost,),:nur)
    )
end

"""
    HQ_const_BG_init()

Field layout for the charm first moment on a prescribed 1+1D background (`matrix1d_visc_HQ_BG!`):
`(alpha, nu^r)`. The background (T, u) is prescribed through `background_fields(tau, r)`, not
evolved. The 2-D counterpart is [`HQ_2p1d_BG_init`](@ref).

`alpha` carries a plain ghost; `nu^r` is odd at the axis.
"""
@inline function HQ_const_BG_init()
    return Fields(
    NDField((:ghost,),(:ghost,),:α), 
    NDField((:odd,),(:ghost,),:nur)
    )
end

"""
    HQ_const_BG_second_moment_init()

[`HQ_const_BG_init`](@ref) plus the charm second moment: `(alpha, nu^r, piQ^l, piQ^phi, Pi_Q)`
(field names `:PiQr, :PiQperp, :PiQ`) on a prescribed background. In radial symmetry the moment is
diagonal in the parallel-transported triad `(l, phihat, etahat)`, so two traceless channels and
the trace `Pi_Q` suffice (`piQ^eta = -(piQ^l + piQ^phi)`). The 2-D counterpart is
[`HQ_2p1d_BG_m2_init`](@ref), which stores four components because the transverse block is not
diagonal.

Layout of `matrix1d_visc_HQ_BG_consistent_m2_5f!`, `matrix1d_visc_HQ_BG_consistent_5f!` and the
not usable `matrix1d_visc_HQ_BG_second_moment!`. Field 5 is the trace `Pi_Q`; the legacy keys
`Meta` (and `Mr`, `Mperp` for fields 3, 4) refer to these fields, not to M^eta_eta.

c_M = 0 by default (`HQ_CM_BACKREACTION`); c_M = D_s/T makes the system elliptic.
"""
@inline function HQ_const_BG_second_moment_init()
    return Fields(
    NDField((:ghost,),(:ghost,),:α), 
    NDField((:odd,),(:ghost,),:nur),
    NDField((:even,),(:ghost,),:PiQr),
    NDField((:even,),(:ghost,),:PiQperp),
    NDField((:even,),(:ghost,),:PiQ)
    )
end

"""
    HQ_2p1d_BG_init()

Field layout for the 2+1D charm sector on a prescribed background: (alpha, nu^x, nu^y).
Cartesian transverse grid, so every boundary is a plain ghost.
"""
@inline function HQ_2p1d_BG_init()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:α),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nuy)
    )
end

"""
    HQ_2p1d_BG_m2_init()

Field layout for the 2+1D charm sector with the consistent second moment
(`matrix2d_HQ_BG_m2!`): (alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy},
piQ^eta_eta, Pi_Q) on a prescribed background.

Every boundary is a plain ghost. Tracelessness is algebraic: `matrix2d_HQ_BG_m2!`
holds `piQ^eta_eta` on its tracelessness value `hq2d_piQ_eta(...)` by a stiff
relaxation in the RHS (rate 100/fm, which limits an explicit integrator's step
to ~0.03 fm).
"""
@inline function HQ_2p1d_BG_m2_init()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:α),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nuy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxx),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQeta),
    NDField((:ghost,:ghost),(:ghost,:ghost),:PiQ)
    )
end

"""
    HQ_viscous_2d()

Field layout for the 10-field 2+1D viscous hydro + passive charm system
(`matrix2d_visc_HQ_BG!`): (T, u^x, u^y, pi^yy, pi^zz, pi^xy, Pi_B, alpha, nu^x, nu^y).
"""
@inline function HQ_viscous_2d()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
    NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piB),
    NDField((:ghost,:ghost),(:ghost,:ghost),:α),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nuy)
    )
end

"""
    HQ_viscous_gamma_1d()

[`HQ_viscous_1d`](@ref) with the charm fugacity carried as `gamma = e^alpha` instead of `alpha`
(`matrix1d_visc_gamma_HQ!`). Equivalent to the alpha form for unit-charge states only; with
charmonia (q = 2) the two differ by ~1 %.
"""
@inline function HQ_viscous_gamma_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:γ),
    NDField((:odd,),(:ghost,),:nur)
    )
end 

"""
    HQ_viscous_no_bulk_1d()

[`HQ_viscous_1d`](@ref) without the bulk pressure: `(T, u^r, pi^phi_phi, pi^eta_eta, alpha, nu^r)`.
Six fields, shear only.
"""
@inline function HQ_viscous_no_bulk_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:α), 
    NDField((:odd,),(:ghost,),:nur)
    )
end 


"""
    HQ_viscous_BG_ideal_current_1d()

Viscous medium plus an IDEAL (non-diffusing) charm fugacity: `(T, u^r, pi^phi_phi, pi^eta_eta, Pi,
alpha)`. No `nu^r` field, so the charm is advected and diluted but does not diffuse (the
`D_s -> 0` limit).
"""
@inline function HQ_viscous_BG_ideal_current_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:α)
    )
end 

"""
    HQ_viscous_BG_ideal_current_density_1d()

[`HQ_viscous_BG_ideal_current_1d`](@ref) with the charm carried as the density `n` instead of the
fugacity `alpha`. The choice of evolved variable changes the time-derivative column of the matrix.
"""
@inline function HQ_viscous_BG_ideal_current_density_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:n)
    )
end 

"""
    HQ_ideal_1d()

The smallest charm layout: ideal medium and an ideal charm fugacity, `(T, u^r, alpha)`. No
dissipative field of any kind.
"""
@inline function HQ_ideal_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:α),
    )
end 

"""
    HQ_ideal_br_p_density_1d()

[`HQ_ideal_1d`](@ref) with the charm as a density: `(T, u^r, n)`.
Its drivers `matrix1d_ideal_HQ_br_p_density!` / `_fugacity!` are not usable.
"""
@inline function HQ_ideal_br_p_density_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:n),
    )
end 

"""
    viscous_1d()

Field layout for the 5-field 1+1D viscous hydro system (`matrix1d_visc!`):
`(T, u^r, pi^phi_phi, pi^eta_eta, Pi)` in cylindrical Milne. The Cartesian sibling is
[`viscous_2d`](@ref).

`pi^phi_phi` and `pi^eta_eta` are MIXED components (`pi^eta_eta = tau^2 pi^{eta eta}`), the same
convention `viscous_2d` uses for `piyy`/`pizz`. On a Bjorken background the Navier-Stokes fixed
points are

    pi^eta_eta -> -4 eta/(3 tau),   pi^phi_phi -> +2 eta/(3 tau),   Pi -> -zeta/tau

The relaxation is Israel-Stewart with `delta_pipi` set by FLUIDUM_DELTA_PIPI (4/3 tau_pi by
default, 0 for plain IS); no `tau_pipi`, no `lambda_piPi`.

`u^r` is odd at the origin, every other field even: r = 0 is a symmetry axis.
"""
@inline function viscous_1d() 
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB)
    )
end 

"""
    viscous_2d()

Field layout for the 7-field 2+1D viscous hydro system (`matrix2d_visc!`):
`(T, u^x, u^y, pi^{yy}, pi^{zz}, pi^{xy}, Pi)` on a transverse-Cartesian, boost-invariant
Milne slice.  The Cartesian sibling of [`viscous_1d`](@ref).

`pi^{yy}` is the Cartesian `pi^y_y` and `pi^{zz}` the MIXED Milne `pi^eta_eta = tau^2 pi^{eta eta}`
-- the same mixed convention `viscous_1d` uses for `piphiphi`/`pietaeta`.

Every boundary is a plain ghost: the transverse grid is Cartesian.
"""
@inline function viscous_2d()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
    NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piB)
    )
end

function fug(T, nhard_profile, r, eos; rdrop = 20, m = 1.5)       
    if r<=rdrop
        return log(nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure))
        else return log(nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure))
    end            
end

function gamma(T, nhard_profile, r, eos; rdrop = 20, m = 1.5)       
    if r<=rdrop
        @show nhard_profile(10)
        @show thermodynamic(T(10),0.0,eos).pressure
        return nhard_profile(r)/(thermodynamic(T(r),0.0,eos).pressure)
        else return nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure)
    end            
end

function gamma_limit_one(T, nhard_profile, r, eos; rdrop = 20, m = 1.5)       
    if r<=rdrop
        return (nhard_profile(r)+r/30*(thermodynamic(T(r),0.0,eos).pressure-nhard_profile(r)))/(thermodynamic(T(r),0.0,eos).pressure)
        else return nhard_profile(rdrop)/(thermodynamic(T(rdrop),0.0,eos).pressure)
    end            
end


function fug_fs(T, ncoll_profile, r, eos;dσ_QQdy,tau_fs, tau0, m = 1.5)       
    return log(density_polar(T, ncoll_profile,r;dσ_QQdy, tau0,tau_fs, m)/(thermodynamic(T,0.0,eos).value))            
end

function fug_fs_rdrop(T, ncoll_profile, r, eos;rdrop = 13, dσ_QQdy, tau_fs, tau0, m = 1.5)       
    if r <= rdrop 
        return fug_fs.(T(r), Ref(ncoll_profile), r, Ref(eos);dσ_QQdy,tau0,tau_fs,m=1.5)   
    else 
        ncoll_rdrop = x -> ncoll_profile(rdrop) 
        return fug_fs(T(rdrop), ncoll_rdrop, rdrop, eos;dσ_QQdy,tau0,tau_fs,m=1.5)
         
    end        
end



function nur_polar(T,ncoll_profile, fonll_profile,r;tau_fs, tau0, dσ_QQdy, m = 1.5) 
    nux_(x,y) = nux(T,ncoll_profile, fonll_profile,x,y; tau_fs, tau0, dσ_QQdy, m) 
    nur_polar = nux_(r,0) 
    return nur_polar
end

"""
Initialize the temperature, the fugacity, and the diffusion current 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, z::TabulatedData{E,F}, cent1::Integer, cent2::Integer; grid_params, initial_params, dσ_QQdy, exp_tail = true, m=1.5) where {A,B,C,D,E,F}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau_fs = initial_params.tau_fs
    tau0 = initial_params.tau0

    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_x = initial_params.norm/tau0, norm_y = 2/tau_fs*dσ_QQdy, exp_tail)
    ccbar_hard = quadgk(x->2*pi*x*tau_fs*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    ccbar_hard_after8 = quadgk(x->2*pi*x*tau_fs*nhard_profile(x),8,rmax,rtol=0.00001)[1]
    Ncoll = quadgk(x->pi*7^2*nhard_profile(x)/(2/tau_fs*dσ_QQdy),0,rmax,rtol=0.00001)[1]
    @show Ncoll
    @show ccbar_hard
    @show ccbar_hard_after8

    eos=Heavy_Quark(readresonancelist(), ccbar_hard)
    fugacity_initial(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop)   
    
    ccbar_thermo,err= quadgk(x->2*pi*x*tau_fs*thermodynamic(temperature_profile(x),fugacity_initial(x),eos.hadron_list).pressure,0,rmax,rtol=0.00001) #integration to obtain the total number of charm
    @show ccbar_thermo
    
    
    fonll_profile = Profiles(z;norm_ = 1e-10, exp_tail = false)
    
    nur_profile(r) = nur_polar.(temperature_profile(r),Ref(nhard_profile), Ref(fonll_profile),r; tau_fs, tau0, dσ_QQdy, m)
    
    fugacity(r) = fug_fs_rdrop(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop, dσ_QQdy, tau_fs, tau0, m = 1.5)  
        
    ccbar_thermo_fs,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=0.00001) #integration to obtain the total number of charm
    
    ccbar_fs,err= quadgk(x->2*pi*x*tau0*density_fs.(temperature_profile(x),Ref(nhard_profile),x,0;dσ_QQdy,tau_fs, tau0, m = 1.5),0,rmax,rtol=0.00001) #integrazione per ottenere il numero totale di charm
    ccbar_fs_after8,err= quadgk(x->2*pi*x*tau0*density_fs.(temperature_profile(x),Ref(nhard_profile),x,0;dσ_QQdy,tau_fs, tau0, m = 1.5),8,rmax,rtol=0.00001) #integrazione per ottenere il numero totale di charm
    
    @show ccbar_thermo_fs
    @show ccbar_fs
    @show ccbar_fs_after8

    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
    set_array!(phi,(x)->nur_profile(x),:nur,disc_fields); #diffusion current initialization
return DiscreteInitialFields(disc,disc_fields,phi), ccbar_thermo_fs
end

"""
Initialize the temperature and the fugacity 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, cent1::Integer, cent2::Integer; grid_params, initial_params, dσ_QQdy, exp_tail = true) where {A,B,C,D}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0

    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_x = initial_params.norm/tau0, norm_y = 2/tau0*dσ_QQdy, exp_tail)
        
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    ccbar_out = quadgk(x->2*pi*x*tau0*nhard_profile(x),8,rmax,rtol=0.00001)[1]
    
    @show ccbar
    @show ccbar_out

    eos=Heavy_Quark(readresonancelist(), ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop)  
    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos.hadron_list).pressure,0,rmax,rtol=0.00001) #integration to obtain the total number of charm
    @show ccbar_thermo

    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
return DiscreteInitialFields(disc,disc_fields,phi), ccbar
end


"""
Get initial temperature and fugacity profiles from tabulated data
"""
function get_initial_T_n_from_tabulated_data(x::TabulatedData{A,B}, y::TabulatedData{C,D}, cent1::Integer, cent2::Integer; grid_params, initial_params,dσ_QQdy, exp_tail = false) where {A,B,C,D}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0
initial_params.rdrop
    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_x = initial_params.norm/tau0,xmax_temp = initial_params.rdrop, xmax_ncoll = initial_params.rdrop, norm_y =2/tau0*dσ_QQdy, exp_tail)

    #@show temperature_profile(8.1)
    #@show nhard_profile(5.1)
    #temperature10 = quadgk(x->temperature_profile(x),0,10,rtol=0.00001)[1]
    #@show temperature10
    #temperaturermax = quadgk(x->temperature_profile(x),0,rmax,rtol=0.00001)[1]
    #@show temperaturermax
    #ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    #@show ccbar
    return Dict{Symbol,Function}(
        :temperature => (x -> temperature_profile(x)),
        :n           => (x -> nhard_profile(x)),
    )
end

"""
Get initial fugacity profile from temperature and density profiles
"""
function get_initial_α_from_T_and_n(temperature_profile, nhard_profile, eos, initial_params)

    fugacity = r -> fug(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop)  
    #fug10 = quadgk(x->fugacity(x),0,10.,rtol=1e-7)[1]
    #@show fugacity(10)
    #@show fug10
    #fugrmax = quadgk(x->fugacity(x),0,30.,rtol=1e-7)[1]
    #@show fugrmax
    return Dict(
        :α => fugacity,
    )
end

function get_initial_γ_from_T_and_n(temperature_profile, nhard_profile, eos, initial_params)

    gamma = r -> gamma_limit_one(temperature_profile, nhard_profile, r, eos.hadron_list;initial_params.rdrop)  
    #fug10 = quadgk(x->fugacity(x),0,10.,rtol=1e-7)[1]
    #@show fugacity(10)
    #@show fug10
    #fugrmax = quadgk(x->fugacity(x),0,30.,rtol=1e-7)[1]
    #@show fugrmax
    return Dict(
        :γ => gamma,
    )
end

"""
initialize fields given an dict of initialization functions
"""
function initialize_fields(init_fun_dict, field_initializer, grid_params)
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints

    oned_visc_hydro = field_initializer
    disc = Fluidum.CartesianDiscretization(OriginInterval(gridpoints, rmax))
    disc_fields = Fluidum.DiscreteFields(oned_visc_hydro(), disc, Float64)
    phi = nothing
    try 
        phi = Fluidum.set_array((x) -> init_fun_dict[:temperature](x), :temperature, disc_fields) #temperature initialization
    catch
        phi = Fluidum.set_array((x) -> init_fun_dict[:α](x), :α, disc_fields) #temperature initialization
    end
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:nu](x), :nur, disc_fields) #diffusion current initialization
    catch 
    end 
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:α](x), :α, disc_fields) #fugacity initialization
    catch
    end 
    try 
        Fluidum.set_array!(phi, x -> init_fun_dict[:γ](x), :γ, disc_fields) # gamma initialization
    catch
    end
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:n](x), :n, disc_fields) #density initialization
    catch
    end
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:PiQr](x), :PiQr, disc_fields)
    catch
        try
            Fluidum.set_array!(phi, x -> init_fun_dict[:Mr](x), :PiQr, disc_fields)
        catch
        end
    end
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:PiQperp](x), :PiQperp, disc_fields)
    catch
        try
            Fluidum.set_array!(phi, x -> init_fun_dict[:Mperp](x), :PiQperp, disc_fields)
        catch
        end
    end
    try
        Fluidum.set_array!(phi, x -> init_fun_dict[:PiQ](x), :PiQ, disc_fields)
    catch
        try
            # legacy key: `:Meta` here means field 5 = the trace Pi_Q, not M^eta_eta; prefer :PiQ.
            Fluidum.set_array!(phi, x -> init_fun_dict[:Meta](x), :PiQ, disc_fields)
        catch
        end
    end
    return Fluidum.DiscreteInitialFields(disc, disc_fields, phi)
end

"""
Initialize fields when only the temperature profile is given
"""
function initialize_fields(x::TabulatedData{A,B}, cent1::Integer, cent2::Integer; grid_params, initial_params, exp_tail = true) where {A,B}
    rmax = grid_params.rmax
    gridpoints = grid_params.gridpoints
    tau0 = initial_params.tau0
    
    temperature_profile = Profiles(x,cent1,cent2; radius = range(0, rmax, gridpoints), norm = initial_params.norm/tau0, exp_tail, temperature_flag = true)         


    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    oned_visc_hydro = Fluidum.viscous_1d()
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    
return DiscreteInitialFields(disc,disc_fields,phi)
end


"""
Use a step-like function for the collision profile if the collision profile is not given from trento
"""
function ncoll_step(r; Ncoll = 1, rdrop = 7, R = 7)       
    if r<=rdrop
        return  Ncoll/(pi*R^2)
        else return 0.01
    end            
end


"""
Initialize fields when only the temperature profile is given, using a step-like function for the collision profile
"""
function initialize_fields(x::TabulatedData{A,B},list; tau0 = 0.4, gridpoints=500,rmax=30,norm_temp=1, norm_coll=1, Ncoll = 1, rdrop = 7) where {A,B}
    temperature_profile = Profiles(x; norm_temp = norm_temp, radius = range(0, rmax, gridpoints)) 
    
    nhard_profile(r) = norm_coll*ncoll_step(r;Ncoll,rdrop)*(temperature_profile(r)/temperature_profile(0.0))^4 
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)

    ccbar_thermo,err= quadgk(x->2*pi*x*tau0*thermodynamic(temperature_profile(x),fugacity(x),eos).pressure,0,rmax,rtol=0.00001) 
    @show ccbar_thermo
    
    n_therm(x) = thermodynamic(temperature_profile(x),fugacity(x),eos).pressure
    
    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
return DiscreteInitialFields(disc,disc_fields,phi)
end

"""
Initialize the temperature, the fugacity, and the diffusion current 
"""
function initialize_fields(x::TabulatedData{A,B}, y::TabulatedData{C,D}, cent1::Integer, cent2::Integer, list; gridpoints=500, rmax=30, norm_temp=1, tau0 = 0.4, dσ_QQdy = 1, exp_tail = true,rdrop = 12,offset = 0.01, m=1.5) where {A,B,C,D}
    temperature_profile, nhard_profile = Profiles(x,y,cent1,cent2; radius = range(0, rmax, gridpoints), norm_temp = norm_temp, norm_coll = 2/tau0*dσ_QQdy, exp_tail, offset=offset)
    
    ccbar = quadgk(x->2*pi*x*tau0*nhard_profile(x),0,rmax,rtol=0.00001)[1]
    @show ccbar

    eos=HadronResonaceGas_ccbar(list, ccbar)
    fugacity(r) = fug(temperature_profile, nhard_profile, r, eos;rdrop)   

    oned_visc_hydro = Fluidum.HQ_viscous_1d()
    disc=CartesianDiscretization(OriginInterval(gridpoints,rmax)) 
    disc_fields = DiscreteFields(oned_visc_hydro,disc,Float64) 
    phi=set_array((x)->temperature_profile(x),:temperature,disc_fields); #temperature initialization
    set_array!(phi,(x)->fugacity(x),:α,disc_fields); #fugacity initialization
    
return DiscreteInitialFields(disc,disc_fields,phi)
end


# ══ Field layouts of the conformal set and the two-way system (closures in Matrix/two_way.jl) ══
# 1+1D: u^r and nu^r odd at the axis, every scalar and every diagonal tensor component even.
# 2+1D: plain ghosts. The S3/S4 1+1D layouts are the existing
# `HQ_const_BG_init` / `HQ_const_BG_second_moment_init`, the S3 2+1D one `HQ_2p1d_BG_init`, the
# generic-EOS S2 ones `viscous_1d` / `viscous_2d`.

"""
    tw_ideal_1d()

1+1D ideal hydro (`tw_ideal_1p1d!`): `(T, u^r)`.
"""
@inline function tw_ideal_1d()
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur)
    )
end

"""
    tw_viscous_conformal_1d()

1+1D viscous hydro with a conformal medium (`tw_viscous_1p1d!` with `ConformalEOS`; no bulk):
`(T, u^r, pi^phi_phi, pi^eta_eta)`, mixed components as in `viscous_1d`.
"""
@inline function tw_viscous_conformal_1d()
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta)
    )
end

"""
    two_way_1d()

The 1+1D two-way system, generic medium EOS (`two_way_1p1d!`), 10 fields:
`(T, u^r, pi^phi_phi, pi^eta_eta, Pi, alpha, nu^r, p_l, p_phi, Pi_Q)`.
Medium as `viscous_1d`; charm as `HQ_const_BG_second_moment_init` (p_l, p_phi on the
parallel-transported triad, p_eta = -(p_l + p_phi)). Charm energy-momentum is in T^{mu nu}_total.
"""
@inline function two_way_1d()
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:piB),
    NDField((:even,),(:ghost,),:α),
    NDField((:odd,),(:ghost,),:nur),
    NDField((:even,),(:ghost,),:PiQr),
    NDField((:even,),(:ghost,),:PiQperp),
    NDField((:even,),(:ghost,),:PiQ)
    )
end

"""
    two_way_conformal_1d()

[`two_way_1d`](@ref) with a conformal medium (no bulk), 9 fields:
`(T, u^r, pi^phi_phi, pi^eta_eta, alpha, nu^r, p_l, p_phi, Pi_Q)`.
"""
@inline function two_way_conformal_1d()
    return Fields(
    NDField((:even,),(:ghost,),:temperature),
    NDField((:odd,),(:ghost,),:ur),
    NDField((:even,),(:ghost,),:piphiphi),
    NDField((:even,),(:ghost,),:pietaeta),
    NDField((:even,),(:ghost,),:α),
    NDField((:odd,),(:ghost,),:nur),
    NDField((:even,),(:ghost,),:PiQr),
    NDField((:even,),(:ghost,),:PiQperp),
    NDField((:even,),(:ghost,),:PiQ)
    )
end

"""
    tw_ideal_2d()

2+1D ideal hydro (`tw_ideal_2p1d!`): `(T, u^x, u^y)`.
"""
@inline function tw_ideal_2d()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
    NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:uy)
    )
end

"""
    tw_viscous_conformal_2d()

2+1D viscous hydro with a conformal medium (no bulk): `(T, u^x, u^y, pi^{yy}, pi^eta_eta, pi^{xy})`
(`pizz` = the mixed pi^eta_eta, as in `viscous_2d`).
"""
@inline function tw_viscous_conformal_2d()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
    NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pixy)
    )
end

"""
    HQ_2p1d_BG_m2_7f_init()

2+1D charm with the consistent second moment on a prescribed background, seven fields
(`tw_hq_m2_bg_2p1d!`): `(alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy}, Pi_Q)`.
piQ^eta_eta is eliminated algebraically (tracelessness + u-orthogonality) instead of being stored
and relaxed as in the 8-field [`HQ_2p1d_BG_m2_init`](@ref).
"""
@inline function HQ_2p1d_BG_m2_7f_init()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:α),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nuy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxx),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:PiQ)
    )
end

"""
    two_way_2d()

The 2+1D two-way system, generic medium EOS (`two_way_2p1d!`), 14 fields:
`(T, u^x, u^y, pi^{yy}, pi^eta_eta, pi^{xy}, Pi, alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy}, Pi_Q)`.
Medium as `viscous_2d`, charm as [`HQ_2p1d_BG_m2_7f_init`](@ref).
"""
@inline function two_way_2d()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
    NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piB),
    NDField((:ghost,:ghost),(:ghost,:ghost),:α),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nuy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxx),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:PiQ)
    )
end

"""
    two_way_conformal_2d()

[`two_way_2d`](@ref) with a conformal medium (no bulk), 13 fields:
`(T, u^x, u^y, pi^{yy}, pi^eta_eta, pi^{xy}, alpha, nu^x, nu^y, piQ^{xx}, piQ^{xy}, piQ^{yy}, Pi_Q)`.
"""
@inline function two_way_conformal_2d()
    return Fields(
    NDField((:ghost,:ghost),(:ghost,:ghost),:temperature),
    NDField((:ghost,:ghost),(:ghost,:ghost),:ux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:uy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:piyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pizz),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pixy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:α),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nux),
    NDField((:ghost,:ghost),(:ghost,:ghost),:nuy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxx),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQxy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:pQyy),
    NDField((:ghost,:ghost),(:ghost,:ghost),:PiQ)
    )
end
