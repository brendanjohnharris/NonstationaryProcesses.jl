module NonstationaryProcesses

using DynamicalSystems
using DifferentialEquations
using Random
using Plots
using Reexport

include("Discontinuous.jl")
include("ParameterProfiles.jl")
include("Process.jl")
include("Processes.jl")
include("Simulations.jl")

end
