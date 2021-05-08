# ------------------------------------------------------------------------------------------------ #
#                                  Sprott's simplest chaotic flow                                  #
# ------------------------------------------------------------------------------------------------ #
# Sprott, or Sprott1997
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

simplestChaoticFlowSim = Process(
    process = simplestChaoticFlow,
    X0 = [0.05, 0.05, 0.05], # As in Sprott's paper
    parameter_profile = ramp,
    parameter_profile_parameters = (2.02, 2.07, 0.0, 10000.0),
    # Parameters should only be between ~2.018 and ~2.082 (otherwise, unbounded)
    # Range can be different when parameters are time varying
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.5,
    tmax = 10000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => false))
export simplestChaoticFlowSim

simplestChaoticFlowArt = Process(
    process = simplestChaoticFlow,
    X0 = [0.05, 0.05, 0.05],
    parameter_profile = ramp,
    parameter_profile_parameters = (2.03, 2.05, 0.0, 5000.0),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.01,
    tmax = 5000.0,
    alg =  RK4(),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-9))
export simplestChaoticFlowArt


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

doublePendulumSim = Process(
    process = doublePendulum,
    X0 = [π, 0.01, 0.0, 0.0], # Two angles and two momenta
    parameter_profile = (ramp, constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = ((0.0, 1.0, 0.0), (1.0,), (1.0,), (1.0,), (1.0,)), # (threshold, baseline, stepHeight)
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.1,
    tmax = 1000.0,
    alg = AutoVern7(Rodas5()))
export doublePendulumSim

cartesianDoublePendulumSim = Process(
    process = cartesianDoublePendulum,
    X0 = [π/2, π/2, 0.0, 0.0], # Two angles and two momenta
    parameter_profile = (ramp, constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = ((-0.015, 1.0, 0.0), (1.0,), (1.0,), (2.0,)), # (threshold, baseline, stepHeight)
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.00001,
    savedt = 0.01,
    tmax = 50.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-20))
export cartesianDoublePendulumSim

cartesianDoublePendulumArt = Process(
    process = cartesianDoublePendulum,
    X0 = [5π/6, π, 0.0, 0.0], # Two angles (from the downward direction) and two momenta
    parameter_profile = (lorentzian, lorentzian, constantParameter, constantParameter),
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 50.0,
    parameter_profile_parameters = ((-1.0, 10.0, 25.0, 1.5), (-1.0, 5.0, 25.0, 2.0), (1.0,), (2.0,)),
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => false))
export cartesianDoublePendulumArt




# ------------------------------------------------------------------------------------------------ #
#                                        Wave-drive Harmonic                                       #
# ------------------------------------------------------------------------------------------------ #
# Sprott, or Chernikov1988
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

waveDrivenHarmonicSim = Process(# parameters:  (ω, ε, k, Ω)
    process = waveDrivenHarmonic,
    X0 = [0.1, 0.1],
    parameter_profile = (constantParameter, constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = ((1π,), (5π,), (1.0π,), (3π,)),
    #parameter_profile_parameters = ((1π,), (10.0,), (1.0π,), (0.9π)),
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.01,
    tmax = 100.0,
    alg = RK4())
export waveDrivenHarmonicSim


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

pulseDrivenHarmonicSim = Process(# parameters:  (ω, ε, k, Ω)
    process = pulseDrivenHarmonic,
    X0 = [0.0, 0.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = ((1π,), (10.0,), (1π,), (1.5π,)),
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 100.0,
    alg = RK4())
export pulseDrivenHarmonicSim



# ------------------------------------------------------------------------------------------------ #
#                              Thomas' cyclically symmetric attractor                              #
# ------------------------------------------------------------------------------------------------ #
# From Sprott, or Elwakil2001

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

thomasCyclicallySymmetricSim = Process(
    process = thomasCyclicallySymmetric,
    X0 = [0.1, 0.0, 0.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (0.18,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 1000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => true))
export thomasCyclicallySymmetricSim

thomasCyclicallySymmetricArt = Process(
    process = thomasCyclicallySymmetric,
    X0 = [0.1, 0.0, 0.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (0.18,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.05,
    tmax = 5000.0,
    alg = RadauIIA3(),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-8))
export thomasCyclicallySymmetricArt


# ------------------------------------------------------------------------------------------------ #
#                                           Double Scroll                                          #
# ------------------------------------------------------------------------------------------------ #
# From Sprott, or Elwakil2001

@inline @inbounds function doubleScroll(X::AbstractArray, 𝑎::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    𝑥̇ = 𝑦
    𝑦̇ = 𝑧
    𝑧̇ = -𝑎(𝑡)*(𝑧 + 𝑦 + 𝑥 - sign(𝑥))
    return SVector{3}(𝑥̇, 𝑦̇, 𝑧̇)
end
function 𝛿(x)
    if x == 0
        return Inf
    else
        return 0
    end
end
@inline @inbounds function doubleScroll_J(X::AbstractArray, 𝑎::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    J = @SMatrix [  0.0   1.0   0.0;
                    0.0   0.0   1.0;
                    -𝑎(𝑡)*(1.0-2.0*𝛿(x))  -𝑎(𝑡)  -𝑎(𝑡)] # Need to worry about integrating delta? Hopefully we never hit 0.0...
end

function doubleScroll(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=doubleScroll_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

doubleScrollSim = Process(
    process = doubleScroll,
    X0 = [0.01, 0.01, 0.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (0.8,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 1000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => true))
export doubleScrollSim

doubleScrollArt = Process(
    process = doubleScroll,
    X0 = [0.01, 0.01, 0.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (0.8,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 5000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => true))
export doubleScrollArt


# ------------------------------------------------------------------------------------------------ #
#                                    Diffusionless Lorenz System                                   #
# ------------------------------------------------------------------------------------------------ #
# Li2014, Sprott2010a (Elegant Chaos), and Schrier2000

@inline @inbounds function diffusionlessLorenz(X::AbstractArray, 𝑎::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    𝑥̇ = 𝑦 - 𝑥
    𝑦̇ = -𝑥*𝑧
    𝑧̇ = 𝑥*𝑦 - 𝑎(𝑡)
    return SVector{3}(𝑥̇, 𝑦̇, 𝑧̇)
end
@inline @inbounds function diffusionlessLorenz_J(X::AbstractArray, 𝑎::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    J = @SMatrix [  -1.0   1.0   0.0;
                    -𝑧     0.0   -𝑥;
                    𝑦      𝑥     0.0]
end

function diffusionlessLorenz(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=diffusionlessLorenz_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

diffusionlessLorenzSim = Process(
    process = diffusionlessLorenz,
    X0 = [1.0, 0.0, 1.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (1.0,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 1000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => false))
export diffusionlessLorenzSim

diffusionlessLorenzArt = Process(
    process = diffusionlessLorenz,
    X0 = [1.0, 0.0, 1.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (1.0,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 5000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => true))
export diffusionlessLorenzArt
