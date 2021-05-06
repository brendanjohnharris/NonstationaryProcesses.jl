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
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=vanderpol_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

# You want to get fancy? Let's store some stuff about the process using multiple dispatch 📜👽
# function vanderpol(s::Symbol)
#     if s == :Equation || s == :equation
#         latexstring("y = \\dot{x}\n\\dot{y} = \\mu(1-x^2)y - x")
#     end
# end




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
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=harmonic_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end




# ------------------------------------------------------------------------------------------------ #
#                                        Harmonic Oscillator                                       #
# ------------------------------------------------------------------------------------------------ #

@inline @inbounds function pitchfork(x::Vector, p::Function, t::Real)
    (μ, α) = p(t)
    dx = μ*x[1] + α*x[1]^3
    return [dx]
end

@inline @inbounds function pitchfork_J(x::Vector, p::Function, t::Real)
    (μ, α) = p(t)
    J = [μ + 3α*x[1]^2]
end

function pitchfork(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=pitchfork_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end




# ------------------------------------------------------------------------------------------------ #
#                                     Skewed Harmonic Oscillator                                   #
# ------------------------------------------------------------------------------------------------ #

# @inline @inbounds function skewedHarmonic(X::AbstractArray, p::Function, t::Real)
#     (ω, κ) = p(t)
#     dX2 = -ω^2.0*X[1] + ω^2*X[1]*exp(-κ*X[1])
#     dX1 = X[2]
#     return SVector{2}(dX1, dX2)
# end

# @inline @inbounds function skewedHarmonic_J(X::AbstractArray, p::Function, t::Real)
#     (ω, κ) = p(t)
#     J = @SMatrix [0.0 1.0; -ω^2 - κ^2*exp(-κ*x) 0.0]
# end

# function skewedHarmonic(P::Process)
#     seed(P.solver_rng)
#     prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=skewedHarmonic_J)
#     sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
# end

