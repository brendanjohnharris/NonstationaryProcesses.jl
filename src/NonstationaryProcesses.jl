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
import Library


function __init__()
    @require BifurcationKit="0f109fa4-8a5d-4b75-95aa-f515264e7665" @eval include("Bifurcations.jl")
    @require DynamicalSystems="61744808-ddfa-5f27-97ff-6e42cc95d634" @eval include("DynamicalSystems.jl")
    @require StatsPlots="f3b207a7-027a-5e70-b257-86293d7955fd" @eval include("Plots/Plotting.jl")
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" @eval include("PyPlotTools.jl")
end

include("Discontinuous.jl")
include("ParameterProfiles.jl")
include("Process.jl")

end
