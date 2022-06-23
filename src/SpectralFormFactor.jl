

#=
function spectral_form_factor(t::Float64, E::Vector{Float64})
    sff_imag = 0.0
    sff_real = 0.0
    for e in E
        sff_imag += sin(e * t * 2*pi)
        sff_real += cos(e * t * 2*pi)
    end

    return sff_real^2 + sff_imag^2
end
=#

#this version is slightly faster
@inline function spectral_form_factor(t::Float64, E::Vector{Float64})
    sff = complex(0.0)
    for e in E
        sff += exp(im * e * t * 2*pi)
    end

    return real(sff)^2 + imag(sff)^2
end


function spectral_form_factor(spect::UnfoldedSpectrum; min = 0.0, max =2.0, grid = 200)
    E =spect.data
    ts = collect(LinRange(min,max, grid))
    sff = zeros(grid) 
    Threads.@threads for i in 1:grid
        sff[i] = spectral_form_factor(ts[i],E)        
    end
    return ts, sff
    
end

SFF = spectral_form_factor