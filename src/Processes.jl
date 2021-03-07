
# A list of process EOM's, Jacobians and simulators (a simulator is a just a function that wraps up parameters for automation)


function seed(theSeed=nothing) # Seed the rng, but return the seed. If no, nothing or NaN argument, randomly seed rng
    if isnothing(theSeed)
        Random.rand(Random.seed!(), UInt64)
    else
        Random.rand(Random.seed!(theSeed), UInt64)
    end
end

# ------------------------------------------------------------------------------------------------ #
#               Define a function which, if it gets a Discontinuity, fills in tstops               #
# ------------------------------------------------------------------------------------------------ #

function dsolve(prob, alg; kwargs...)
    if typeof(prob.p) <: Discontinuous
        DifferentialEquations.solve(prob, alg; kwargs..., tstops=sort(collect(prob.p.d))) # May need to check tstops isn't in args in the future
    else
        DifferentialEquations.solve(prob, alg; kwargs...)
    end
end

function tuplef2ftuple(f, params)
    if typeof(f) <: Tuple
        ps = Vector{Function}(undef, length(f))
        for i = 1:length(f)
            ps[i] = f[i](params[i]...)
        end
        p(t) = map((x, g) -> g(x), fill(t, length(ps)), ps) # Something like that
    else
        p = f(params...)
    end
end


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

function vanderpol(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=vanderpol_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
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

function henon(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=henon_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end



# ------------------------------------------------------------------------------------------------ #
#                                        Harmonic Oscillator                                       #
# ------------------------------------------------------------------------------------------------ #

@inline @inbounds function harmonic(X::AbstractArray, ω::Function, t::Real)
    dX2 = -ω(t)^2.0*X[1]
    dX1 = X[2]
    return SVector{2}(dX1, dX2)
end

@inline @inbounds function harmonic_J(X::AbstractArray, ω::Function, t::Real)
    J = @SMatrix [0.0 1.0; -ω(t)^2 0.0]
end

function harmonic(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=harmonic_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end


# ------------------------------------------------------------------------------------------------ #
#                                            Noisy Sine                                            #
# ------------------------------------------------------------------------------------------------ #
function noisySine(P::Process)
    seed(P.solver_rng)
    sol = [sin(t + asin(P.X0...)) + P.parameter_profile(P.parameter_profile_parameters...)(t)*randn() for t in P.t0:P.savedt:P.tmax]
end

