using .PyCall
using .DifferentialEquations
const cn = PyNULL()
const signal = PyNULL()
run(`$(PyCall.python) -m pip install colorednoise`)
copy!(cn, pyimport("colorednoise"))
copy!(signal, pyimport_conda("scipy.signal", "scipy"))

function seed(theSeed=nothing) # Seed the rng, but return the seed. If no, nothing or NaN argument, randomly seed rng
    if isnothing(theSeed)
        theSeed = abs(Random.rand(Int64))
    end
    Random.seed!(theSeed)
    return theSeed
end

# ------------------------------------------------------------------------------------------------ #
#               Define a function which, if it gets a Discontinuity, fills in tstops               #
# ------------------------------------------------------------------------------------------------ #

function dsolve(prob, alg; kwargs...)
    if prob.p isa Discontinuous
        DifferentialEquations.solve(prob, alg; kwargs..., tstops=sort(collect(prob.p.d))) # May need to check tstops isn't in args in the future
    else
        DifferentialEquations.solve(prob, alg; kwargs...)
    end
end

"""
Whip up an ODEproblem from a Process
"""
function process2problem(P::Process)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))
end
export process2problem

function tuplef2ftuple(f, params)
    # turn a tuple of functions into a function of tuples
    if all(isempty.(params)) # The f's are just functions on their own, no need to add parameters
        # Be warned that you can't mix these; either use all parameter functions, or all standard functions. Don't be greedy.
        if f isa Tuple
            pp(t) = [x(t) for x in f]
        else
            pp = f
        end
        return pp
    end
    if f isa Tuple
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
include("Processes/Transforms.jl")
