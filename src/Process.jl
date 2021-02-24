# How best to structure the I/O format and system definitions
# The DynamicalSystems package already has its own neat method of storing system information and time series, but not the simulation parameters. You want to be able to save this to a file that can be read ealsewhere, so we don't want any Julia-specific types in our dict (probably don't want this regardless, for data persistence).
# I think we'll want a simple dictionary, with enough fields that you can reconstruct the function calls that generate a time series dataset.

Base.@kwdef mutable struct Process
    process::Union{Function, Nothing} = nothing
    X0::Union{AbstractArray, Nothing} = nothing
    parameter_function::Union{String, Nothing, Tuple} = nothing # Can be a tuple of strings, if the system has more than one parameter
    parameter_function_parameters::Union{Tuple, Nothing} = nothing # Can be a tuple of tuples
    # parameter_offset::Union{Real, Nothing, Tuple} = nothing
    # parameter_scale::Union{Real, Nothing, Tuple} = nothing
    transient::Union{Real, Nothing} = nothing
    tmax::Union{Real, Nothing} = nothing
    dt::Union{Real, Nothing} = nothing
    t0::Union{Real, Nothing} = nothing
    parameter_rng::Union{Real, Nothing} = nothing
    solver_rng::Union{Real, Nothing} = nothing
    solver::Union{String, Nothing} = nothing
    relative_tolerance::Union{Real, Nothing} = nothing
    inventory_id::Int = abs(rand(Int, 1)[1]) # Just a unique number for this simulation
    solver_opts::Union{Tuple, Nothing} = nothing
    # Add any other solver options as they become necessary
end


# ------------------------------------------------------------------------------------------------ #
#               A function to handle simulations that are specified with a Process type            #
# ------------------------------------------------------------------------------------------------ #
function simulate(P::Process)
    # Generate the parameter function, and save a parameter vector
    P.parameter_rng = copy(seed(P.parameter_rng))
    parameter_function = eval.(Meta.parse.(P.parameter_function))
    if isnothing(P.parameter_function_parameters)
        parameter_function_parameters = fill([], length(parameter_function))
    else
        parameter_function_parameters = P.parameter_function_parameters
    end
    if typeof(parameter_function) <: Tuple
        ps = Vector{Function}(undef, length(parameter_function))
        for i = 1:length(parameter_function)
            ps[i] = parameter_function[i](parameter_function_parameters[i]...)
        end
        p(t) = map((x, g) -> g(x), fill(t, length(ps)), ps) # Something like that
    else
        p = parameter_function(parameter_function_parameters...)
    end
    T = P.t0:P.dt:P.tmax
    seed(P.parameter_rng) # So the function and the evaluated parameters have the same rng state, if there is anything stochastic in the function itself
    parameters = p(T)
    if isnothing(P.solver_opts)
        solver_opts = ()
    else
        solver_opts = P.solver_opts
    end
    data, metadata = P.process(X0=P.X0,
                                p=p,
                                T= (P.t0, P.t0, P.dt, P.tmax),
                                solver=eval(Meta.parse(P.solver)),
                                reltol=P.relative_tolerance,
                                rngseed=P.solver_rng, solver_opts...)

    # If any of the fields of P are nothing, we check to see if the simulation gives us any answers in metadata
    fns = fieldnames(typeof(P))
    for fn = 1:length(fns)
        subFn = fns[fn] # Throws a syntax error otherwise
        if isnothing(getproperty(P, subFn))
            setproperty!(P, subFn, getproperty(metadata, subFn))
        end
    end

    return (Dict(:Trajectory => data, :Parameters => parameters), P)
end
export simulate

function timeseries(s::Dict, dim::Real=1)
    s[:Trajectory][:, dim]
end
function timeseries(s::Dict, dim::Union{Vector, Tuple})
    map(x -> timeseries(s, x), dim)
end
function timeseries(s::Tuple, dim::Union{Real, Vector, Tuple}=1)
    timeseries(s[1], dim) # You gave the metadata as well
end
export timeseries

function Discrete() end
export Discrete # Hotfix
