module NonstationaryProcesses

using Requires
using Random
using Reexport
using StatsBase
using Tullio
using FFTW
@reexport using NonstationaryProcessesBase
using Distributions

function __init__()
    @require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" @eval begin
        include("../ext/DifferentialEquationsExt/DifferentialEquationsExt.jl")
    end
end

include("ARMA.jl")
include("Noise.jl")
include("Signals.jl")

end # module
