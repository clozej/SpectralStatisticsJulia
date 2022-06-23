
include("Probability.jl")

function level_spacing(spect::UnfoldedSpectrum; n::Int = 1)
    e = spect.data
    s = e[1 + n : end] .- e[1 : end-n]
    return s
end

# make a macro that adds pdf to statitic name and computes it!
function level_spacing_pdf(spect::UnfoldedSpectrum, bins::Vector{T}; n::Int = 1) where T
    s = level_spacing(spect; n=n)
    return pdf_hist(s, bins)
end

function level_spacing_cdf(spect::UnfoldedSpectrum, pts::Vector{T}; n::Int = 1) where T
    s = level_spacing(spect; n=n)
    return cdf(s, bins)
end

function level_spacing_u_cdf(spect::UnfoldedSpectrum, pts::Vector{T}; n::Int = 1) where T
    s = level_spacing(spect; n=n)
    return u_cdf(s, bins)
end


function level_spacing_ratio(spect::UnfoldedSpectrum, shift::Int=1, n::Int = 1)
    s = level_spacing(spect, n=n)
    shifted = circshift(s,-shift)
    r = s[1:end-shift] ./ shifted[1:end-shift]
    return r
end


