using LaTeXStrings
# A list of process EOM's, Jacobians and simulators (a simulator is a just a function that wraps up parameters for automation)


function seed(theSeed=nothing) # Seed the rng, but return the seed. If no, nothing or NaN argument, randomly seed rng
    if isnothing(theSeed)
        Random.rand(Random.seed!(), UInt64)
    else
        Random.rand(Random.seed!(theSeed), UInt64)
    end
end

# ------------------------------------------------------------------------------------------------ #
#               Define a function which, if it gets a Discontinuity, fills in tstops               #
# ------------------------------------------------------------------------------------------------ #

function dsolve(prob, alg; kwargs...)
    if typeof(prob.p) <: Discontinuous
        DifferentialEquations.solve(prob, alg; kwargs..., tstops=sort(collect(prob.p.d))) # May need to check tstops isn't in args in the future
    else
        DifferentialEquations.solve(prob, alg; kwargs...)
    end
end

function tuplef2ftuple(f, params)
    # turn a tuple of functions into a function of tuples
    if typeof(f) <: Tuple
        ps = Vector{Function}(undef, length(f))
        for i = 1:length(f)
            ps[i] = f[i](params[i]...)
        end
        p(t) = map((x, g) -> g(x), fill(t, length(ps)), ps) # Something like that
    else
        p = f(params...)
    end
end
export tuplef2ftuple

include("Processes/ARMA.jl")
include("Processes/ChaoticFlows.jl")
include("Processes/ChaoticMaps.jl")
include("Processes/DeterministicFlows.jl")
include("Processes/Noise.jl")
include("Processes/Signals.jl")
include("Processes/Stochastic.jl")
