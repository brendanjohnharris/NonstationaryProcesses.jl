module NonstationaryProcesses

using Random
using Reexport
using StatsBase
using Tullio
using FFTW
using StaticArrays
@reexport using NonstationaryProcessesBase
using Distributions

# ODE/SDE-based simulation machinery. Formerly a Requires/weakdep extension; now a
# direct dependency. OrdinaryDiffEq 7 (the SciMLBase-3 build) re-exports SciMLBase
# (ODEProblem/solve) plus the Verner solvers (AutoVern, Vern), but no longer
# bundles RK4 (LowOrderRK), RadauIIA3 (FIRK), FunctionMap, or Rodas5 (Rosenbrock),
# so those solver packages are pulled in explicitly. StochasticDiffEq provides the SDE solver (EM) and its `solve` method.
using OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK   # RK4
using OrdinaryDiffEqFIRK         # RadauIIA3
using OrdinaryDiffEqFunctionMap  # FunctionMap
using OrdinaryDiffEqRosenbrock   # Rodas5
using StochasticDiffEq           # EM, SDEProblem
import NonstationaryProcessesBase: dsolve, process2problem, process2solution, odeproblem

"""Solve, filling in tstops when the parameters carry Discontinuities."""
function dsolve(prob, alg; kwargs...)
    if prob.p isa Discontinuous
        solve(prob, alg; kwargs..., tstops=sort(collect(prob.p.d)))
    else
        solve(prob, alg; kwargs...)
    end
end

"""Whip up a problem from a Process."""
function process2problem(P::Process, jac=nothing)
    args = [P.process, P.X0, (P.transient_t0, P.tmax),
        tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters)]
    isnothing(jac) || append!(args, jac)
    if getalg(P) isa FunctionMap
        DiscreteProblem(args...)
    else
        odeproblem(args...)
    end
end

process2solution(P::Process, jac=nothing) =
    dsolve(process2problem(P, jac), P.alg; dt=P.dt, saveat=P.savedt, P.solver_opts...)

function odeproblem(f, X0, ts, ps; jac=nothing)
    odefunction = ODEFunction{true}(f; jac)
    return ODEProblem(odefunction, X0, ts, ps)
end

include("ARMA.jl")
include("Noise.jl")
include("Signals.jl")
include("ChaoticFlows.jl")
include("ChaoticMaps.jl")
include("DeterministicFlows.jl")
include("Stochastic.jl")
include("Transforms.jl")

end # module
