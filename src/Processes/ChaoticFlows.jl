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

"""
    Sprott's simplest chaotic flow

Sprott, or Sprott1997
"""
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

simplestChaoticFlowVis = Process(
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
export simplestChaoticFlowVis


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

"""
    Your classic double pendulum, with Hamiltonian equations of motion

    https://scienceworld.wolfram.com/physics/DoublePendulum.html
"""
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
    parameter_profile = (constantParameter, constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = ((1.0,), (1.0,), (1.0,), (1.0,)), # (ℓ₁, ℓ₂, m₁, m₂)
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.00001,
    savedt = 0.1,
    tmax = 100.0,
    alg = AutoVern9(Rodas5()),
    solver_opts =  Dict(:adaptive => true, :reltol => 1e-12, :abstol => 1e-12))
export cartesianDoublePendulumSim

cartesianDoublePendulumVis = Process(
    process = cartesianDoublePendulum,
    X0 = [5π/6, π, 0.0, 0.0], # Two angles (from the downward direction) and two momenta
    parameter_profile = (lorentzian, lorentzian, constantParameter, constantParameter),
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 50.0,
    parameter_profile_parameters = ((-1.0, 10.0, 25.0, 1.5), (-1.0, 5.0, 25.0, 2.0), (1.0,), (2.0,)),
    alg = AutoVern9(Rodas5()),
    solver_opts = Dict(:adaptive => false))
export cartesianDoublePendulumVis


cartesianDoublePendulumVis2 = Process(
    process = cartesianDoublePendulum,
    X0 = [5π/6, π, 0.0, 0.0], # Two angles (from the downward direction) and two momenta
    parameter_profile = (lorentzian, lorentzian, constantParameter, constantParameter),
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 100.0,
    parameter_profile_parameters = ((-1.0, 18.0, 50.0, 1.5), (-1.0, 9.0, 50.0, 2.0), (1.0,), (2.0,)),
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => false))
export cartesianDoublePendulumVis2



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

"""
    An harmonic oscillator driven by a plane wave propagating in the x-direction

Sprott, or Chernikov1988
"""
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


# @inline @inbounds function pulseDrivenHarmonic(X::AbstractArray, p::Function, t::Real)
#     (ω, ε, k, Ω) = p(t)
#     dX2 = -ω^2.0*X[1] + ε*sin(k*X[1] - Ω*t)^2
#     dX1 = X[2]
#     return SVector{2}(dX1, dX2)
# end

# @inline @inbounds function pulseDrivenHarmonic_J(X::AbstractArray, p::Function, t::Real)
#     (ω, ε, k, Ω) = p(t)
#     J = @SMatrix [0.0 1.0; 2*ε*k*cos(k*X[1] - Ω*t)-ω^2 0.0]
# end


# """
#     An harmonic oscillator driven by directional 'pulses'
# """
# function pulseDrivenHarmonic(P::Process)
#     seed(P.solver_rng)
#     prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))#, jac=waveDrivenHarmonic_J)
#     sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
# end

# pulseDrivenHarmonicSim = Process(# parameters:  (ω, ε, k, Ω)
#     process = pulseDrivenHarmonic,
#     X0 = [0.0, 0.0],
#     parameter_profile = (constantParameter, constantParameter, constantParameter, constantParameter),
#     parameter_profile_parameters = ((1π,), (10.0,), (1π,), (1.5π,)),
#     transient_t0 = -10.0,
#     t0 = 0.0,
#     dt = 0.001,
#     savedt = 0.01,
#     tmax = 100.0,
#     alg = RK4())
# export pulseDrivenHarmonicSim



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

"""
    Thomas' cyclically symmetric attractor

From Sprott, or Thomas1999
"""
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

thomasCyclicallySymmetricVis = Process(
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
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))
export thomasCyclicallySymmetricVis



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

"""
    Double scroll attractor

From Sprott, or Elwakil2001
"""
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

doubleScrollVis = Process(
    process = doubleScroll,
    X0 = [0.01, 0.01, 0.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (0.8,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 5000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))
export doubleScrollVis



@inline @inbounds function lorenz(dX, X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    (𝜎, 𝑟, 𝑏) = p(𝑡)
    dX[1] = 𝜎*(𝑦 - 𝑥)
    dX[2] = -𝑥*𝑧 + 𝑟*𝑥 - 𝑦
    dX[3] = 𝑥*𝑦 - 𝑏*𝑧
end

@inline @inbounds function lorenz_J(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    (𝜎, 𝑟, 𝑏) = p(𝑡)
    J = @SMatrix [  -𝜎      𝜎       0.0;
                    𝑟-𝑧     -1.0    -𝑥;
                    𝑦       𝑥       -𝑏]
end

"""
    Lorenz system

Sprott, Shin1998, Balcerzak2018
"""
function lorenz(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=lorenz_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end


lorenzSim = Process(
    process = lorenz,
    X0 = [0.0, -0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3), # Sprott's recomendation
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 500.0,
    alg = AutoVern9(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-10, :abstol => 1e-10, :maxiters=>1e7))
export lorenzSim

lorenzVis= Process(
    process = lorenz,
    X0 = [0.0, 0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.0005,
    savedt = 0.001,
    tmax = 1000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-15))
export lorenzVis



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

"""
    Diffusionless Lorenz system

See Li2014, Sprott2010a (Elegant Chaos), and Schrier2000
"""
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

diffusionlessLorenzVis = Process(
    process = diffusionlessLorenz,
    X0 = [1.0, 0.0, 1.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (1.0,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.025,
    tmax = 20000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-8))
export diffusionlessLorenzVis



@inline @inbounds function piecewiseLinearHyperchaos(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧, 𝑢) = X
    (𝑎, 𝑏) = p(𝑡)

    𝑥̇ = 𝑦 - 𝑥
    𝑦̇ = -𝑧*sign(𝑥) + 𝑢
    𝑧̇ = abs(𝑥) - 𝑎
    𝑢̇ = -𝑏*𝑦

    return SVector{4}(𝑥̇, 𝑦̇, 𝑧̇, 𝑢̇)
end
@inline @inbounds function piecewiseLinearHyperchaos_J(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧, 𝑢) = X
    (𝑎, 𝑏) = p(𝑡)
    J = @SMatrix [  -1.0        1.0        0.0          0.0;
                    -2*𝛿(𝑧)     0.0        -sign(𝑥)     1.0; # ? Again, hope we miss the zero
                    sign(𝑥)     0.0        0.0          0.0;
                    0.0         -𝑏         0.0          0.0]
end

"""
    A simple, piecewise-linear, hyperchaotic flow

Li2014
"""
function piecewiseLinearHyperchaos(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=piecewiseLinearHyperchaos_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

piecewiseLinearHyperchaosSim = Process(
    process = piecewiseLinearHyperchaos,
    X0 = [0.2, 0.0, 1.0, 0.0],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = (1.0, 0.25), # * 𝑎 is an amplitude, so set to 1 without loss of generality
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 1000.0,
    alg = Vern9(),
    solver_opts = Dict(:adaptive => false, :reltol => 1e-10, :abstol => 1e-10))
export piecewiseLinearHyperchaosSim

piecewiseLinearHyperchaosVis = Process(
    process = piecewiseLinearHyperchaos,
    X0 = [0.01, 0.01, 0.01, 0.01],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = (1.0, 0.25),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 15000.0,
    alg = AutoVern9(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-14))
export piecewiseLinearHyperchaosVis




@inline @inbounds function simplifiedLorenz4D(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧, 𝑢) = X
    (𝑎, 𝑏) = p(𝑡)

    𝑥̇ = 𝑦 - 𝑥
    𝑦̇ = -𝑥*𝑧 + 𝑢
    𝑧̇ = 𝑥*𝑦 - 𝑎
    𝑢̇ = -𝑏*𝑦

    return SVector{4}(𝑥̇, 𝑦̇, 𝑧̇, 𝑢̇)
end
@inline @inbounds function simplifiedLorenz4D_J(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧, 𝑢) = X
    (𝑎, 𝑏) = p(𝑡)
    J = @SMatrix [  -1.0        1.0        0.0          0.0;
                    -𝑧          0.0        -𝑥           1.0;
                    𝑦           𝑥          0.0          0.0;
                    0.0         -𝑏         0.0          0.0]
end

"""
    The simplified 4D Lorenz system, whose attractors are exclusively hidden and can be hyperchaotic

Li2014a and Gao2006
"""
function simplifiedLorenz4D(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=simplifiedLorenz4D_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

simplifiedLorenz4DSim = Process(
    process = simplifiedLorenz4D,
    X0 = [2.0, 4.0, 0.0, 0.0],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = (2.6, 0.44),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.05,
    tmax = 1000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => false))
export simplifiedLorenz4DSim

simplifiedLorenz4DVis = Process(
    process = simplifiedLorenz4D,
    X0 = [2.0, 4.0, 0.0, 0.0],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = (2.6, 0.45),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 7500.0,
    alg = Vern9(),
    solver_opts = Dict(:adaptive => false))#, :reltol => 1e-12))
export simplifiedLorenz4DVis



@inline @inbounds function chensSystem(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    (𝑎, 𝑏, 𝑐) = p(𝑡)

    𝑥̇ = 𝑎*(𝑦 - 𝑥)
    𝑦̇ = (𝑐 - 𝑎)*𝑥 - 𝑥*𝑧 + 𝑐*𝑦
    𝑧̇ = 𝑥*𝑦 - 𝑏*𝑧

    return SVector{3}(𝑥̇, 𝑦̇, 𝑧̇)
end
@inline @inbounds function chensSystem_J(X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    (𝑎, 𝑏, 𝑐) = p(𝑡)
    J = @SMatrix [-𝑎        𝑎       0.0;
                  𝑐-𝑎-𝑧     𝑐       -𝑥;
                  𝑦         𝑥       -𝑏]
end


"""
    Yet another chaotic attractor

Chen1999, Song2004
"""
function chensSystem(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=chensSystem_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end
export chensSystem
yetAnotherChaoticAttractor = chensSystem
export yetAnotherChaoticAttractor

chensSystemSim = Process(
    process = chensSystem,
    X0 = [-3.0, 2.0, 20.0],
    parameter_profile = (ramp, ramp, ramp),
    parameter_profile_parameters = ((42.0, 46.0, 0.0, 500.0), (4.0, 11.0, 0.0, 500.0), (28.0, 38.0, 0.0, 500.0)),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.025,
    tmax = 500.0,
    alg = AutoVern9(Rodas5()),
    solver_opts = Dict(:adaptive => true, :abs_tol=>1e-10, :rel_tol=>1e-10))
export chensSystemSim

chensSystemVis = Process(
    process = chensSystem,
    X0 = [-3.0, 2.0, 20.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (35.0, 8.0/3.0, 28.0),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 1000.0,
    alg = Vern9(),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-9))
export chensSystemVis
