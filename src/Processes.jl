using LaTeXStrings
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
export tuplef2ftuple


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
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=henon_J)
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
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=harmonic_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end


# ------------------------------------------------------------------------------------------------ #
#                                            Noisy Sine                                            #
# ------------------------------------------------------------------------------------------------ #
function noisySine(P::Process)
    seed(P.solver_rng)
    sol = [sin(t + asin(P.X0...)) + P.parameter_profile(P.parameter_profile_parameters...)(t)*randn() for t in P.transient_t0:P.savedt:P.tmax]
end



# ------------------------------------------------------------------------------------------------ #
#                                      Noisy Trendy Scaly Sine                                     #
# ------------------------------------------------------------------------------------------------ #
function noisyShiftyScalySine(P::Process)
    # Parameters (η, C, A)
    seed(P.solver_rng)
    (η, C, A) = [P.parameter_profile[x](P.parameter_profile_parameters[x]...) for x in 1:length(P.parameter_profile)]
    sol = [A(t)*sin(t + asin(P.X0...)) + η(t)*randn() + C(t) for t in P.transient_t0:P.savedt:P.tmax]
end

function shcalySine(P::Process)
    seed(P.solver_rng)
    A = P.parameter_profile(P.parameter_profile_parameters...)
    sol = [A(t)*(sin(t + asin(P.X0...)) + 1.0*randn() + 1.0) for t in P.transient_t0:P.savedt:P.tmax]
end



# ------------------------------------------------------------------------------------------------ #
#                                          Double Pendulum                                         #
# ------------------------------------------------------------------------------------------------ #

polarReduce(x::Real) = x - 2π*(x÷π);
polarReduce(X::AbstractArray) = X .- 2π.*(X.÷π); export polarReduce

@inline @inbounds function doublePendulum(X::AbstractArray, p, t::Real)
    (θ₁, θ₂, 𝑝₁, 𝑝₂) = X
    (ℓ₁, ℓ₂, m₁, m₂) = p(t)
    𝑔 = 1.0

    C₁(θ₁, θ₂, 𝑝₁, 𝑝₂) = (𝑝₁*𝑝₂*sin(θ₁ - θ₂))/(ℓ₁*ℓ₂*(m₁ + m₂*sin(θ₁ - θ₂)^2))
    C₂(θ₁, θ₂, 𝑝₁, 𝑝₂) = ((ℓ₂^2*m₂*𝑝₁^2 + ℓ₁^2*(m₁ + m₂)*𝑝₂^2 - ℓ₁*ℓ₂*m₂*𝑝₁*𝑝₂*cos(θ₁ - θ₂))/(2*ℓ₁^2*ℓ₂^2*(m₁ + m₂*sin(θ₁ - θ₂)^2)^2))*(sin(2*(θ₁ - θ₂)))

    dθ₁ = (ℓ₂*𝑝₁ - ℓ₁*𝑝₂*cos(θ₁ - θ₂))/(ℓ₁^2*ℓ₂*(m₁ + m₂*sin(θ₁ - θ₂)^2))
    dθ₂ = (ℓ₁*(m₁ + m₂)*𝑝₂ - ℓ₂*m₂*𝑝₁*cos(θ₁ - θ₂))/(ℓ₁*ℓ₂^2*m₂*(m₁ + m₂*sin(θ₁ - θ₂)^2))

    d𝑝₁ = -(m₁ + m₂)*𝑔*ℓ₁*sin(θ₁) - C₁(θ₁, θ₂, 𝑝₁, 𝑝₂) + C₂(θ₁, θ₂, 𝑝₁, 𝑝₂)
    d𝑝₂ = -m₂*𝑔*ℓ₂*sin(θ₂) + C₁(θ₁, θ₂, 𝑝₁, 𝑝₂) - C₂(θ₁, θ₂, 𝑝₁, 𝑝₂)

    return SVector{4}(dθ₁, dθ₂, d𝑝₁, d𝑝₂)
end
function doublePendulum(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

function cartesianDoublePendulum(P::Process)
    # Should the initial conditions be set in cartesian coordinates?
    sol = doublePendulum(P(process=doublePendulum))
    p = parameters(P, transient=true)
    x₁ = p[1, :].*sin.(sol[1, :])
    y₁ = .-p[1, :].*cos.(sol[1, :])
    x₂ = x₁ .+ p[2, :].*sin.(sol[2, :])
    y₂ = y₁ .- p[2, :].*cos.(sol[2, :])
    sol = hcat(x₁, y₁, x₂, y₂)
end


# ------------------------------------------------------------------------------------------------ #
#                                        Wave-drive Harmonic                                       #
# ------------------------------------------------------------------------------------------------ #
# https://doi.org/10.1063/1.881159
@inline @inbounds function waveDrivenHarmonic(X::AbstractArray, p::Function, t::Real)
    (ω, ε, k, Ω) = p(t)
    dX2 = -ω^2.0*X[1] + ε*sin(k*X[1] - Ω*t)
    dX1 = X[2]
    return SVector{2}(dX1, dX2)
end

@inline @inbounds function waveDrivenHarmonic_J(X::AbstractArray, p::Function, t::Real)
    (ω, ε, k, Ω) = p(t)
    J = @SMatrix [0.0 1.0; ε*k*cos(k*X[1] - Ω*t)-ω^2 0.0]
end

function waveDrivenHarmonic(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=waveDrivenHarmonic_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

# ------------------------------------------------------------------------------------------------ #
#                                        Pulse-drive Harmonic                                      #
# ------------------------------------------------------------------------------------------------ #
@inline @inbounds function pulseDrivenHarmonic(X::AbstractArray, p::Function, t::Real)
    (ω, ε, k, Ω) = p(t)
    dX2 = -ω^2.0*X[1] + ε*sin(k*X[1] - Ω*t)^2
    dX1 = X[2]
    return SVector{2}(dX1, dX2)
end

# @inline @inbounds function pulseDrivenHarmonic_J(X::AbstractArray, p::Function, t::Real)
#     (ω, ε, k, Ω) = p(t)
#     J = @SMatrix [0.0 1.0; 2*ε*k*cos(k*X[1] - Ω*t)-ω^2 0.0]
# end

function pulseDrivenHarmonic(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))#, jac=waveDrivenHarmonic_J)
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



# ------------------------------------------------------------------------------------------------ #
#                                   Skewed Pitchfork Bifurcation                                   #
# ------------------------------------------------------------------------------------------------ #

@inline @inbounds function skewedGaussianQuadratic(x::AbstractArray, p::Function, t::Real)
    (μ, η, κ) = p(t)
    x = μ*x[1] + μ*x[1]*exp(-κ*x[1])*(κ*x[1]-2.0)
    return SVector{1}(x)
end

@inline @inbounds function skewedGaussianQuadratic_σ(σ::AbstractArray, p::Function, t::Real)
    (μ, η, κ) = p(t)
    σ = η
end

function skewedGaussianQuadratic(P::Process)
    seed(P.solver_rng)
    prob = SDEProblem(P.process, skewedGaussianQuadratic_σ, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

