module DifferentialEquationsExt
using ..DifferentialEquations
using StatsBase
using Tullio
using ..DifferentialEquations.DiffEqBase
using ..DifferentialEquations.OrdinaryDiffEq
using StaticArrays
using FFTW
using ..NonstationaryProcesses
import NonstationaryProcessesBase: dsolve, process2problem, process2solution, odeproblem


include("./ChaoticFlows.jl")
include("./ChaoticMaps.jl")
include("./DeterministicFlows.jl")
include("./Stochastic.jl")
include("./Transforms.jl")
end # module
