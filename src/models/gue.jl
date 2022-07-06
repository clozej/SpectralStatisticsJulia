using SpecialFunctions 
using Base.MathConstants


struct GUE <: Model end


#first argument of all functions is the variable


function level_spacing_pdf(gue::GUE, s; n=1)
    #Wigner surmises for higher orders Wen-Jia Rao
    beta = 2.0 #level repulsion
    a = 0.5*n*(n+1.0)*beta + n -1
    A =  (gamma(0.5*a + 1.0)/gamma(0.5*a + 0.5) )^2.0 
    C = (2.0*(gamma(0.5*a + 1.0))^(a + 1.0)) / ((gamma(0.5*a + 0.5))^(a + 2.0))
    return @. C * s^a * exp(-A * s^2.0)
end

function level_spacing_cdf(gue::GUE, s; n=1)
    beta = 2.0 #level repulsion
    a = 0.5*n*(n+1.0)*beta + n -1
    A =  (gamma(0.5*a + 1.0)/gamma(0.5*a + 0.5) )^2.0 
    C = (2.0*(gamma(0.5*a + 1.0))^(a + 1.0)) / ((gamma(0.5*a + 0.5))^(a + 2.0))
    return @. 0.5*C*s^(1.0 + a)*(A*s^2.0)^(0.5*(-1.0 - a))*(gamma((1.0 + a)*0.5) - gamma((1.0 + a)*0.5, A*s^2.0))

end

function level_spacing_u(gue::GUE, s)
    return @. (2.0 / pi) * arccos(sqrt(1.0 - level_spacing_cdf(gue, s)))    
end

function number_variance(gue::GUE, l)
    arg = 2.0  * pi .* l
    si, ci = sinint.(arg), cosint.(arg) 
    nv_gue = @. 1.0/pi^2.0 * (log(arg) + eulergamma + 1.0 - cos(arg) - ci) + l * (1.0 - 2.0/pi * si)    
    return nv_gue
end 

function rigidity(gue::GUE, l)
    beta = 2.0
    nv = @. 2.0/(beta*pi^2.0)*(log(beta*pi*l) + eulergamma + 1.0)
    #nv = @. 1.0/(pi^2.0)*(log(beta*pi*l) + eulergamma + 1.0)#approximate formula O(1/(2*pi*l)) 
    return @. 0.5 * nv - 9.0/(4.0*beta*pi^2.0)
end


function spectral_form_factor(gue::GUE, t)
    sff = @. abs.(t) 
    idx = t.> 1.0
    sff[idx] .= 1.0
    return sff
end