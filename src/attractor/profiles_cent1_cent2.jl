# The (cent1, cent2) signatures that phd-git calls, next to upstream's (cent) ones.

"""
Initialize a profile when both the entropy and the density of binary collision profiles are given"""
function Profiles(x::TabulatedData{A,B}, y::TabulatedData{A,B}, cent1::Integer, cent2::Integer; radius, norm_x, norm_y, exp_tail = true,offset_temp=0.01, xmax_temp = 8, xmax_ncoll = 5, offset_ncoll = 0.) where {A,B}
    #entropy profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm_x)
    
    imax = findlast(≤(xmax_temp), r)   # last index where r[i] ≤ xmax_temp
    r = isnothing(imax) ? r[0:0] : r[1:imax]

    entropy_profile = entropy_profile[1:length(r)]
    #@show offset_temp 
    #@show xmax_temp
    #@show offset_ncoll
    #@show xmax_ncoll

    temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile)
    temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
    #temperature_funct = Dierckx.Spline1D(r, temperature_profile, k=1)
    temp_exp_funct = nothing
    if exp_tail == true
        temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = xmax_temp, offset = offset_temp)
        # temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax = 8, offset = 0.01)
        temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat())
    end
    #ncoll profile
    r, ncoll_profile = get_profile(y, cent1, cent2; norm = norm_y)  

    imax = findlast(≤(xmax_temp), r)   # last index where r[i] ≤ xmax_temp
    r = isnothing(imax) ? r[0:0] : r[1:imax]

    ncoll_profile = ncoll_profile[1:length(r)]

    ncoll_funct = linear_interpolation(r, ncoll_profile; extrapolation_bc=Flat())

    ncoll_exp_funct = nothing
    if exp_tail == true
        ncoll_exp = exponential_tail_pointlike.(Ref(ncoll_funct), radius; xmax = xmax_ncoll, offset= offset_ncoll)
        ncoll_exp_funct = linear_interpolation(radius, ncoll_exp; extrapolation_bc=Flat())
    end
    
    if exp_tail == true
        return temp_exp_funct, ncoll_exp_funct
    else
        return temperature_funct, ncoll_funct
    end
end

"""
Initialize when only one profile is given"""
function Profiles(x::TabulatedData{A,B}, cent1::Integer, cent2::Integer; radius, norm, xmax = 8, exp_tail = true, temperature_flag = true, offset = 0.01) where {A,B}
    #entopy profile
    r, entropy_profile = get_profile(x, cent1, cent2; norm = norm)
    
    if temperature_flag == true
        temperature_profile = InverseFunction(x->pressure_derivative(x,Val(1),FluiduMEoS())).(entropy_profile)  
        temperature_funct = linear_interpolation(r, temperature_profile; extrapolation_bc=Flat())
    else
        temperature_funct = linear_interpolation(r, entropy_profile; extrapolation_bc=Flat())
    end

    if exp_tail == true
        temp_exp = exponential_tail_pointlike.(Ref(temperature_funct), radius; xmax, offset)
        temp_exp_funct = linear_interpolation(radius, temp_exp; extrapolation_bc=Flat())
        return temp_exp_funct
    else
        return temperature_funct
    end
end
