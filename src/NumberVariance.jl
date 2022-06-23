include("DataSample.jl")
using StatsBase

function number_variance(L::Float64, E)
    Ave1 = 0.0
    Ave2 = 0.0
    Ave3 = 0.0
    Ave4 = 0.0
    j = 1
    k = 1
    x = E[1]
    largest_energy = E[end - Int(ceil(L)+100)]
    while x < largest_energy
        #move index k
        while E[k] < x+L
            k = k+1    
        end

        d1 = E[j] - x
        d2 = E[k] - (x+L)
        cn = k - j
        if d1 < d2
            x = E[j]
            s = d1
            j = j + 1
        else
            x = E[k] - L
            s = d2
            k = k + 1
        end
        Ave1 = Ave1 + s*cn
        Ave2 = Ave2 + s*cn^2
    end
    s = largest_energy - E[1]
    Ave1 = Ave1/s
    Ave2 = Ave2/s
    AveSig = (Ave2 - Ave1^2)
    return AveSig
end

function number_variance(spect::UnfoldedSpectrum, x::Vector{T}) where T 
    E = spect.data
    Ls = x
    sz = length(x)
    nvs = zeros(sz) 
    Threads.@threads for i in 1:sz
        nvs[i] = number_variance(Ls[i],E)        
    end
    return nvs
end


NV = number_variance