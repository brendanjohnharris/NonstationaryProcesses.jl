module NonstationaryProcesses

using Reexport
@reexport using DifferentialEquations
using Requires
using Random
using StaticArrays

include("Discontinuous.jl")
include("ParameterProfiles.jl")
include("Process.jl")
include("Processes.jl")
include("Plotting.jl")

using BifurcationKit
using Setfield
include("Bifurcations.jl")
# function __init__()
#     @require BifurcationKit="0f109fa4-8a5d-4b75-95aa-f515264e7665" @eval include("Bifurcations.jl")
# end

end
