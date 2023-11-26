using StatsBase

@inline @inbounds function skewedGaussianQuadratic(x::AbstractArray, p::Function, t::Real)
    (μ, η, κ) = p(t)
    x = μ * x[1] + μ * x[1] * exp(-κ * x[1]) * (κ * x[1] - 2.0)
    return SVector{1}(x)
end

@inline @inbounds function skewedGaussianQuadratic_σ(σ::AbstractArray, p::Function, t::Real)
    (μ, η, κ) = p(t)
    σ = η
end

"""Skewed Pitchfork Bifurcation"""
function skewedGaussianQuadratic(P::Process)
    seed(P.solver_rng)
    prob = SDEProblem(P.process, skewedGaussianQuadratic_σ, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))
    sol = dsolve(prob, P.alg; dt=P.dt, saveat=P.savedt, P.solver_opts...)
end

skewedGaussianQuadraticSim = Process(
    process=skewedGaussianQuadratic,
    X0=[0.0],
    parameter_profile=(constantParameter, constantParameter, rampInterval),
    parameter_profile_parameters=((-1.0,), (2.0,), (0.01, 0.05, 0.0, 1000.0)),
    transient_t0=-100.0,
    t0=0.0,
    dt=0.001,
    savedt=0.001,
    tmax=1000.0,
    alg=EM())
export skewedGaussianQuadraticSim


"""A system combining multiple classes of features. Combine Chen's system with gaussian noise and an exponential measurement function"""
function seeIfItSticks(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(chensSystem, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile[1:3], P.parameter_profile_parameters[1:3]), jac=chensSystem_J)
    sol = dsolve(prob, P.alg; dt=P.dt, saveat=P.savedt, P.solver_opts...)
    η = parameterseries(P, p=[4], transient=true)
    α = parameterseries(P, p=[5], transient=true)
    sol = sol .+ η .* randn(size(sol))
    sol = (sol .* exp.(α .* (sol .- mean(sol, dims=2))))'
end

seeIfItSticksSim = Process(
    process=seeIfItSticks,
    X0=[-3.0, 2.0, 20.0],
    parameter_profile=(ramp, ramp, ramp, constant, constant),
    parameter_profile_parameters=((42.0, 42.0, 0.0, 500.0), (11.0, 11.0, 0.0, 500.0), (28.0, 28.0, 0.0, 500.0), 2.0, 0.01), # (𝑎, 𝑏, 𝑐, η, α)
    transient_t0=-100.0,
    t0=0.0,
    dt=0.001,
    savedt=0.025,
    tmax=500.0,
    alg=RK4(),
    solver_opts=Dict(:adaptive => false))
export seeIfItSticksSim
