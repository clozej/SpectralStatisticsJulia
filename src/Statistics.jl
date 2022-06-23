include("AbstractTypes.jl")
include("Probability.jl")
include("LevelSpacings.jl")
include("NumberVariance.jl")

#abstract type Statistic end

function compute_statistic(data::DataSample, statistic::F, args... ; kwargs...) where F
    stat = Symbol(statistic)
    expr = :($stat) 
    return eval(expr)(data, args... ; kwargs...)
end


function compute_statistic(model::M, statistic::F, args... ; kwargs...) where {M, F}
    stat = Symbol(statistic)
    f = getfield(model, stat)
    return f(args... ; kwargs...)
end
