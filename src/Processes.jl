
# A list of process EOM's, Jacobians and simulators (a simulator is a just a function that wraps up parameters for automation)


function seed(theSeed=nothing) # Seed the rng, but return the seed. If no or NaN argument, randomly seed rng
    if isnothing(theSeed)
        Random.rand(Random.seed!(), UInt64)
    else
        Random.rand(Random.seed!(theSeed), UInt64)
    end
end

# Hijack DynamicalSystems trajectory so that (ONLY) if it gets a Discontinuity, it fills in tstops
function DynamicalSystems.trajectory(ds, T; args...)
    if typeof(ds.p) <: Discontinuous
        DynamicalSystems.trajectory(ds, T; tstops=sort(collect(ds.p.d)), args...) # May need to check tstops isn't in args in the future
    else
        DynamicalSystems.trajectory(ds, T; args...)
    end
end

#abstract type ProcessSimulator <: Function end


# ------------------------------------------------------------------------------------------------ #
#                                      Van der Pol Oscillator                                      #
# ------------------------------------------------------------------------------------------------ #
@inline @inbounds function vanderpol(X::AbstractArray, μ::Function, t::Real)
    dX1 = X[2]
    dX2 = μ(t).*(1-X[1].^2).*X[2] - X[1]
    return SVector{2}(dX1, dX2)
end

@inline @inbounds function vanderpol_J(X::AbstractArray, μ::Function, t::Real)
    J = @SMatrix [  0                  1;
                    -2*μ(t)*X[1]*X[2]-1      μ(t)*(1-X[1]^2)]
end

function vanderpol(; X0::AbstractArray,
                     p::Function,
                     T::NTuple{4,Real},
                     solver,
                     reltol::Real,
                     rngseed::Union{UInt, Nothing}, solver_opts...)

    # Do the actual simulation
    rngseed = seed(rngseed)
    ds = ContinuousDynamicalSystem(vanderpol, X0, p, vanderpol_J, t0=T[2])
    data = trajectory(ds, T[4]; dt=T[3], Ttr=T[1], reltol=reltol, alg=solver, solver_opts...)

    # Save the results
    metadata = Process(
    process = vanderpol,
    X0 = X0,
    transient = T[1],
    t0 = T[2],
    dt = T[3],
    tmax = T[4],
    solver_rng = rngseed,
    relative_tolerance = reltol)

    return (data, metadata)
end

# ------------------------------------------------------------------------------------------------ #
#                                             Henon map                                            #
# ------------------------------------------------------------------------------------------------ #
@inline @inbounds function henon(x, p, n) # In this case, p is a tuple of functions
    dx1 = 1.0 - p(n)[1]*x[1]^2 + x[2]
    dx2 = p(n)[2]*x[1]
    return SVector{2}(dx1, dx2)
end
# Jacobian:
@inline @inbounds function henon_J(x, p, n)
    J = @SMatrix [-2*p(n)[1]*x[1]      1.0;
                p(n)[2]            0.0]
    return J
end


function henon(; X0::AbstractArray,
    p::Function,
    T::NTuple{4,Real},
    solver,
    reltol::Real,
    rngseed::Union{UInt, Nothing}, solver_opts...)

    # Do the actual simulation
    rngseed = seed(rngseed)
    ds = DiscreteDynamicalSystem(henon, X0, p, henon_J, t0=T[2])
    show(ds)
    data = trajectory(ds, T[4]; dt=T[3], Ttr=T[1])

    # Save the results
    metadata = Process(
    process = henon,
    X0 = X0,
    transient = T[1],
    t0 = T[2],
    dt = T[3],
    tmax = T[4],
    solver_rng = rngseed,
    relative_tolerance = reltol)

    return (data, metadata)
end
