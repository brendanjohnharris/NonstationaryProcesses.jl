module NonstationaryProcesses

using Reexport
@reexport using DifferentialEquations
using Requires
using Random
using StaticArrays
using PyCall

include("Discontinuous.jl")
include("ParameterProfiles.jl")
include("Process.jl")
include("Processes.jl")
include("Plotting.jl")

using BifurcationKit
using Setfield
include("Bifurcations.jl")

const cn = PyNULL()
const signal = PyNULL()
function __init__()
    #@require BifurcationKit="0f109fa4-8a5d-4b75-95aa-f515264e7665" @eval include("Bifurcations.jl")
    @require DynamicalSystems="61744808-ddfa-5f27-97ff-6e42cc95d634" @eval include("DynamicalSystems.jl")

    run(`$(PyCall.python) -m pip install colorednoise`)
    copy!(cn, pyimport("colorednoise"))
    copy!(signal, pyimport_conda("scipy.signal", "scipy"))
end

end
