using DifferentialEquations
using DifferentialEquations.OrdinaryDiffEq

"""Define a function that, if it gets a Discontinuity, fills in tstops"""
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
function process2problem(P::Process, jac=nothing)
    args = [P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters)]
    isnothing(jac) || append!(ags, jac)
    if getalg(P) isa FunctionMap
        DiscreteProblem(args...)
    else
        ODEProblem(args...)
    end
end
export process2problem

process2solution(P::Process, jac=nothing) = dsolve(process2problem(P, jac), P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
export process2solution

include("Processes/ChaoticFlows.jl")
include("Processes/ChaoticMaps.jl")
include("Processes/DeterministicFlows.jl")
include("Processes/Stochastic.jl")
include("Processes/Transforms.jl")
