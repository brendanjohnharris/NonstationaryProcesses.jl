module NonstationaryProcesses

using Requires
using Random
using StaticArrays
using PyCall
using Distributions
using Tullio
using FFTW
using Setfield
using StatsBase


function __init__()
    @require BifurcationKit="0f109fa4-8a5d-4b75-95aa-f515264e7665" @eval include("Bifurcations.jl")
    @require DynamicalSystems="61744808-ddfa-5f27-97ff-6e42cc95d634" @eval include("DynamicalSystems.jl")
    @require DifferentialEquations="0c46a032-eb83-5123-abaf-570d42b7fbaa" @eval include("DifferentialEquations.jl")
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" @eval include("PyPlotTools.jl")
end

include("Discontinuous.jl")
include("ParameterProfiles.jl")
include("Process.jl")
include("Plotting/Plotting.jl")
include("Processes.jl")

end
