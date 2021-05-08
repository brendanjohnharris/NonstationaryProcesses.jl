using DimensionalData
using DelimitedFiles
using Dates
# How best to structure the I/O format and system definitions...

Base.@kwdef mutable struct Process # Not ensemble
    process = nothing
    parameter_profile::Union{Function, Tuple, Array} = constantParameter # Can be a tuple of symbols, if the system has more than one parameter
    parameter_profile_parameters::Union{Tuple, Array} = [0] # Can be a tuple of tuples
    X0::Vector = [nothing]
    t0::Union{Float64, Int64} = 0.0
    transient_t0::Union{Float64, Int64} = t0 # t0 will always take precedence
    dt::Union{Float64, Int64} = 0.001
    savedt::Union{Float64, Int64} = 0.01
    tmax::Union{Float64, Int64} = 100.0
    alg::Union{SciMLBase.SciMLAlgorithm, Nothing} = nothing
    solver_opts::Dict = Dict(:adaptive=>false)
    #parameter_rng::UInt64 = seed()
    solver_rng::Int64 = seed()
    id::Int64 = abs(rand(Int64)) # Just a unique number for this simulation
    date::String = string(Dates.now())
    solution = nothing
end
function Process(D::Dict)
    for s ∈ setdiff(keys(D), (:date, :solution))
        if typeof(D[s]) <: String
            D[s] = eval(Meta.parse(D[s]))
        end
    end
    Process(;D...)
end
function (P::Process)(;kwargs...)
    # Can use field aliases here
    kwargs = Dict(kwargs)
    repalias!(kwargs, process_aliases)

    P2 = deepcopy(P)
    [setfield!(P2, x, y) for (x, y) in kwargs]
    setfield!(P2, :solution, nothing) # You've changed some parameters, so the solution is no longer valid
    setfield!(P2, :id, abs(rand(Int64))) # New id, yeah?
    setfield!(P2, :solver_rng, seed()) # New random seed, yeah?
    setfield!(P2, :date, string(Dates.now())) # New datetime, yeah?
    return P2
end
export Process

# ------------------------------------------------------------------------------------------------ #
#                                           Field Aliases                                          #
# ------------------------------------------------------------------------------------------------ #
process_aliases = Dict(
    :process =>                     [:sim, :system, :processes],
    :parameter_profile =>           [:profile, :profiles, :𝑝, :𝑃],
    :parameter_profile_parameters =>[:parameters, :ps, :params, :param, :parameter,
                                     :profile_parameters, :parameterprofileparameters,
                                     :profileparameters, :𝔓, :𝔭],
    :X0 =>                          [:initial_conditions, :X, :X_0, :X₀, :𝑥₀, :𝑋₀, :𝑥0,
                                     :𝑋0],
    :transient_t0 =>                [:transient, :cutoff, :tₜ, :𝑡ₜ, :tt],
    :t0 =>                          [:tstart, :t₀, :𝑡₀],
    :dt =>                          [:δt, :𝛿t, :δ𝑡, :𝛿𝑡],
    :savedt =>                      [:save_dt, :save_Δt, :save_δt, :save_𝛥t, :save_𝛿t, :save_Δ𝑡,
                                     :save_δ𝑡, :save_𝛥𝑡, :save_𝛿𝑡, :Δt, :𝛥t, :Δ𝑡, :Δ𝑡],
    :tmax =>                        [:t_max, :T, :Tmax, :T_max, :𝑇],
    :alg =>                         [:algorithm, :solver],
    :solver_opts =>                 [:opts, :solopts, :sol_opts, :solveropts],
    :solver_rng =>                  [:rng, :rngseed, :rng_seed, :solverrng],
    :id =>                          [:identifier, :inventory_id],
    :date =>                        [:time, :datetime],
    :solution =>                    [:sol, :result, :output]
)
function repalias!(D, aliai::Dict)
    for d ∈ keys(D)
        for a ∈ keys(aliai)
            if d ∈ aliai[a]
                D[a] = pop!(D, d)
            end
        end
    end
end

# ------------------------------------------------------------------------------------------------ #
#               A function to handle simulations that are specified with a Process type            #
# ------------------------------------------------------------------------------------------------ #
function solution!(P::Process) # vars::Tuple=Tuple(1:size(P.X0)[1])
    if isnothing(P.solution)
        @info "Solving for the $(getprocess(P)) process ($(getid(P)))"
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
function timeseries!(P::Process, dim=1:length(getX0(P)); transient::Bool=false)
    x = timeseries(solution!(P), dim)
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

function parameter_functions(P::Process)
    if all(isempty.(getps(P)))
        getprofiles(P)
    else
        if typeof(P.parameter_profile) <: Union{Tuple, Vector}
            [P.parameter_profile[x](P.parameter_profile_parameters[x]...) for x in 1:length(P.parameter_profile)]
        else
            P.parameter_profile(P.parameter_profile_parameters...)
        end
    end
end
export parameter_functions

function parameterseries(P::Process; p=nothing, kwargs...)
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
export parameterseries

# Access the fields of a process with functions
for field ∈ keys(process_aliases)
    f = Symbol(:get, field)
    eval(quote
        $f(P::Process) = P.$field; export $f
    end)
    for field_alias ∈ process_aliases[field]
        fa = Symbol(:get, field_alias)
        eval(quote
            $fa = $f; #export fa
        end)
    end
end

function forcevec(x)
    if !(typeof(x) <: Union{AbstractArray, Tuple})
        x = [x]
    else
        x
    end
end
export forcevec

function forcemat(x)
    if typeof(x) <: AbstractVector
        x = reshape(x, :, 1)
    end
    return x
end
export forcemat

function trimtransient(P::Process)
    if !isempty(P.solution)
        P.solution = timeseries(P, transient=false)
    end
    P.transient_t0 = P.t0
    return P
end

function saveTimeseries!(P::Process, folder::String="./", delim::Char=','; transient::Bool=true, fileroot="timeseries")
    X = timeseries(P, transient=transient)
    if !transient
        P = trimtransient(P)
    end
    filename = joinpath(folder, fileroot*"_"*string(getid(P))*".csv")
    P.solution = abspath(folder)
    #P.solution = nothing
    @info "Saving time-series data to $filename"
    writedlm(filename, X, delim)
end
export saveTimeseries!

function timeseries(P::Process, dim=1:length(getX0(P)); folder::Union{String, Bool}=(typeof(getsolution(P)) <: String), kwargs...)
    if typeof(folder) <: Bool && folder
        if typeof(getsolution(P)) <: String
            folder = getsolution(P)
        else
            folder = "./"
        end
    end
    if typeof(folder) <: String
        filename = filter(x->occursin("timeseries_"*string(getid(P)), x), readdir(folder))
    else
        return timeseries!(P, dim; kwargs...)
    end
    if isempty(filename)
        return timeseries!(P, dim; kwargs...)
    end
    filename = filename[1]
    @info "Loading time-series data from $filename"
    filename = joinpath(folder, filename)
    P.solution = readdlm(filename, ',', Float64)
    return timeseries!(P, dim; kwargs...)
end
