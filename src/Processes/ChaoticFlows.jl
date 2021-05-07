# ------------------------------------------------------------------------------------------------ #
#                                  Sprott's simplest chaotic flow                                  #
# ------------------------------------------------------------------------------------------------ #
# https://doi.org/10.1016/S0375-9601(97)00088-1
@inline @inbounds function simplestChaoticFlow(X::AbstractArray, 𝐴::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    𝑥̇ = 𝑦
    𝑦̇ = 𝑧
    𝑧̇ = -𝐴(𝑡)*𝑧 + 𝑦^2 - 𝑥
    return SVector{3}(𝑥̇, 𝑦̇, 𝑧̇)
end

@inline @inbounds function simplestChaoticFlow_J(X::AbstractArray, 𝐴::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    J = @SMatrix [  1.0     0.0     0.0;
                    0.0     1.0     0.0;
                    2.0*𝑦   -𝐴(𝑡)   -1.0]
end

function simplestChaoticFlow(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=simplestChaoticFlow_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end




# ------------------------------------------------------------------------------------------------ #
#                                          Double Pendulum                                         #
# ------------------------------------------------------------------------------------------------ #
# https://scienceworld.wolfram.com/physics/DoublePendulum.html

polarReduce(x::Real) = x - 2π*(x÷π);
polarReduce(X::AbstractArray) = X .- 2π.*(X.÷π);

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
    p = parameterseries(P, transient=true)
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
#                              Thomas' cyclically symmetric attractor                              #
# ------------------------------------------------------------------------------------------------ #
# From Sprott, or https://doi.org/10.1142/S0218127499001383

@inline @inbounds function thomasCyclicallySymmetric(X::AbstractArray, 𝑏::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    𝑥̇ = -𝑏(𝑡)*𝑥 + sin(𝑦)
    𝑦̇ = -𝑏(𝑡)*𝑦 + sin(𝑧)
    𝑧̇ = -𝑏(𝑡)*𝑧 + sin(𝑥)
    return SVector{3}(𝑥̇, 𝑦̇, 𝑧̇)
end

@inline @inbounds function thomasCyclicallySymmetric_J(X::AbstractArray, 𝑏::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    J = @SMatrix [  -𝑏(𝑡)   cos(𝑦)   0.0;
                    0.0     -𝑏(𝑡)    cos(𝑧);
                    cos(𝑥)  0.0      -𝑏(𝑡)]
end

function thomasCyclicallySymmetric(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=thomasCyclicallySymmetric_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end
