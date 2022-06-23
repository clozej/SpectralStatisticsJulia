

@inline function length_spectrum(l::Float64, rho_fluct::Vector{Float64}, ks::Vector{Float64})
    l_spec = complex(0.0)
    for (rho, k)  in zip(rho_fluct, ks) 
        l_spec += rho * exp(im * k * l)
    end
    norm = length(rho_fluct)
    return abs(l_spec)/norm
end


function length_spectrum(rho_fluct::Vector{Float64}, ks::Vector{Float64}; min = 0.0, max =5.0, grid = 200)
    ls = collect(LinRange(min,max, grid))
    l_spec = zeros(grid) 
    Threads.@threads for i in 1:grid
        l_spec[i] = length_spectrum(ls[i],rho_fluct,ks)        
    end
    return ls, l_spec
end