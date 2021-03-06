module NonstationaryProcesses

using Reexport
#@reexport using DynamicalSystems
@reexport using DifferentialEquations
@reexport using Plots
using Random
using StaticArrays

include("Discontinuous.jl")
include("ParameterProfiles.jl")
include("Process.jl")
include("Processes.jl")
include("Simulations.jl")

end
