using SpecialFunctions 
using Base.MathConstants


struct Poisson <: Model end


#first argument of all functions is the variable


function level_spacing_pdf(poisson::Poisson, s ; n=1)
    return @. n^n/factorial(n-1) * s^(n-1) * exp(-n*s)
end

function level_spacing_cdf(poisson::Poisson, s; n=1)
    return @. 1.0 - gamma(n, n*s)/gamma(n)
end

function level_spacing_u(poisson::Poisson, s)
    return @. (2.0 / pi) * arccos(sqrt(1.0 - level_spacing_cdf(poisson, s)))    
end

function number_variance(poisson::Poisson, l) 
  return l
end 

function rigidity(poisson::Poisson, l) 
    return l ./ 15.0
  end 

function spectral_form_factor(poisson::Poisson, t)
    return ones(length(t))
end