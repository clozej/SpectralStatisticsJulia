include("AbstractTypes.jl")

using SpecialFunctions 
using Base.MathConstants


struct GOE_statistics <: Model
    level_spacing_pdf
    level_spacing_cdf
    level_spacing_u
    number_variance
    spectral_form_factor

    function GOE_statistics()
        return new(P, W, U, num_var, SFF)
    end
end


#first argument of all functions is the variable


function P(s)
    #Wigner-Dyson
    a = gamma((3.0) / (2.0))^2.0
    return @. 2.0 * a * s * exp(-a * s^2.0)
end

function W(s)
    a = gamma((3.0) / (2.0))^2.0
    return @. 1.0 - exp(-a * s^2.0)
end

function U(s)
    return @. (2.0 / pi) * arccos(sqrt(1.0 - W(s)))    
end

function num_var(l) 
    pl = pi .* l
    pl2 = 2.0*pi .* l 
    si1, ci1 = sinint.(pl), cosint.(pl) 
    si2, ci2 = sinint.(pl2), cosint.(pl2) 
    
    return @. 2.0/pi^2.0*(log(pl2) + eulergamma + 1.0 +0.5*si1^2 - pi/2.0 * si1 - cos(pl2) - ci2 + pi^2.0 * l * (1.0-2.0/pi * si2))

end 

function SFF(t)
    sff = @. 2.0 * t - t * log(1.0 + 2.0 * t)
    idx = t.> 1.0
    sff[idx] = @. 2.0 - t[idx] * log((1.0 + 2.0 * t[idx])/(2.0 * t[idx] - 1.0))    
    return sff
end