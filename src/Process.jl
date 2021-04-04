using DimensionalData
# How best to structure the I/O format and system definitions...

Base.@kwdef mutable struct Process # Not ensemble
    process = nothing
    parameter_profile::Union{Function, Tuple, Array} = constantParameter # Can be a tuple of symbols, if the system has more than one parameter
    parameter_profile_parameters::Union{Tuple, Array} = [0] # Can be a tuple of tuples
    X0::Vector = [nothing]
    transient_t0::Union{Float64, Int64} = t0 # t0 will always take precedence
    t0::Union{Float64, Int64} = -10.0
    dt::Union{Float64, Int64} = 0.001
    savedt::Union{Float64, Int64} = 0.01
    tmax::Union{Float64, Int64} = 100.0
    alg::Union{SciMLBase.SciMLAlgorithm, Nothing} = nothing
    solver_opts::Dict = Dict(:adaptive=>false)
    #parameter_rng::UInt64 = seed()
    solver_rng::UInt64 = seed()
    inventory_id::UInt64 = abs(rand(UInt64, 1)[1]) # Just a unique number for this simulation
    solution = nothing
end
function (P::Process)(;kwargs...) # Cleaner way to do this with constructors???
    P2 = deepcopy(P)
    [setfield!(P2, x, y) for (x, y) in kwargs]
    setfield!(P2, :solution, nothing) # You've changed some parameters, so the solution is no longer valid
    setfield!(P2, :solver_rng, seed()) # New random seed, yeah?
    return P2
end
export Process

# ------------------------------------------------------------------------------------------------ #
#               A function to handle simulations that are specified with a Process type            #
# ------------------------------------------------------------------------------------------------ #
function solution!(P::Process) # vars::Tuple=Tuple(1:size(P.X0)[1])
    if isnothing(P.solution)
        P.solution = P.process(P)
    end
    x = P.solution
end
export solution
function simulate(P::Process)
    P2 = deepcopy(P)
    P2.solution = solution!(P2)
    return P2
end
simulate!(P::Process) = begin P.solution = solution!(P); return nothing end
export simulate
export simulate!

# --------- Might have various solution types, so timeseries gets all of them as an array -------- #
timeseries(s::SciMLBase.AbstractTimeseriesSolution, dim::Real) = s[dim, :]
timeseries(s::SciMLBase.AbstractTimeseriesSolution, dim::Union{Vector, UnitRange}=1:size(s.u[1], 1)) = copy(s[dim, :]')
function timeseries(s::AbstractArray, dim::Union{Vector, UnitRange, Real}=1:size(s, 2))
    if typeof(s) <: Vector
        if length(dim) != 1 || dim[1] != 1
            error("Cannot index the second dimension of the input, which is a vector")
        end
        s[:]
    else
        s[:, dim]
    end
end
function timeseries(P::Process, args...; transient::Bool=false)
    x = timeseries(solution!(P), args...)
    if transient
        idxs = 1:length(times(P, transient=true))
    else
        idxs = (length(P.transient_t0:P.savedt:P.t0)):1:length(times(P, transient=true))
    end
    saveTimes = (P.transient_t0:P.savedt:P.tmax)[idxs]
    if size(x, 2) > 1
        x = DimArray(x[idxs, :], (Ti(saveTimes), Dim{:Variable}(1:size(x, 2))))
    else
        x = DimArray(x[idxs], (Ti(saveTimes),))
    end
end
# function timeseries(s::Tuple, dim::Union{Real, Vector, Tuple}=1)
#     timeseries(s[1], dim) # You gave the metadata as well
# end
export timeseries

# ------------------------------------------------------------------------------------------------ #
#                              Convenient access to some solution info                             #
# ------------------------------------------------------------------------------------------------ #
function times(P::Process; transient::Bool=false)
    if transient
        P.transient_t0:P.savedt:P.tmax
    else
        P.t0:P.savedt:P.tmax
    end
end
export times

parameter_function(P::Process) = tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters)
export parameter_function

function parameters(P::Process; p=nothing, kwargs...)
    ps = hcat(parameter_function(P).(times(P; kwargs...))...)
    if size(ps, 1) == 1 # This 1 × N array, which should be a vector
        ps = ps[:]
    end
    if isnothing(p) || (typeof(ps)<:Vector)
        return ps
    else
        return ps[p, :]
    end
end
export parameters
