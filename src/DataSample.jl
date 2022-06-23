include("AbstractTypes.jl")

struct RealSpectrum <: DataSample 
    data::Vector{Float64}
end

struct UnfoldedSpectrum <: DataSample 
    data::Vector{Float64}
end

struct WavenumberSpectrum <: DataSample 
    data::Vector{Float64}
end

struct ComplexSpectrum <: DataSample 
    data::Vector{ComplexF64}
end