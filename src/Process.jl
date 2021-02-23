# How best to structure the I/O format and system definitions
# The DynamicalSystems package already has its own neat method of storing system information and time series, but not the simulation parameters. You want to be able to save this to a file that can be read ealsewhere, so we don't want any Julia-specific types in our dict (probably don't want this regardless, for data persistence).
# I think we'll want a simple dictionary, with enough fields that you can reconstruct the function calls that generate a time series dataset.

Base.@kwdef mutable struct Process
    process::Union{Function, Nothing} = nothing
    X0::Union{AbstractArray, Nothing} = nothing
    parameter_function::Union{String, Nothing, Tuple{String}} = nothing # Can be a tuple of strings, if the system has more than one parameter
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
    # Add any other solver options as they become necessary
end


# ------------------------------------------------------------------------------------------------ #
#               A function to handle simulations that are specified with a Process type            #
# ------------------------------------------------------------------------------------------------ #
function simulate(P::Process)
    # Generate the parameter function, and save a parameter vector
    P.parameter_rng = copy(seed(P.parameter_rng))
    parameter_function = eval(Meta.parse(P.parameter_function))
    if typeof(parameter_function) <: Tuple
        for i = 1:length(parameter_function)
            ps(i) = parameter_function(i)(P.parameter_function_parameters(i)...)
        end
        p(t) = [ps...](t) # Something like that
    else
        p = parameter_function(P.parameter_function_parameters...)
    end
    T = P.t0:P.dt:P.tmax
    seed(P.parameter_rng) # So the function and the evaluated parameters have the same rng state, if there is anything stochastic in the function itself
    parameters = p(T)

    data, metadata = P.process(X0=P.X0,
                                p=p,
                                T= (P.t0, P.t0, P.dt, P.tmax),
                                solver=eval(Meta.parse(P.solver)),
                                reltol=P.relative_tolerance,
                                rngseed=P.solver_rng)

    # If any of the fields of P are nothing, we check to see if the simulation gives us any answers in metadata
    fns = fieldnames(typeof(P))
    for fn = 1:length(fns)
        subFn = fns[fn] # Throws a syntax error otherwise
        if isnothing(getproperty(P, subFn))
            setproperty!(P, subFn, getproperty(metadata, subFn))
        end
    end

    return (Dict(:Trajectory => data, :Parameters => parameters), metadata)
end
export simulate