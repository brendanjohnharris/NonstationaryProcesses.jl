@inline @inbounds function forcedvanderpol(X::AbstractArray, p::Function, t::Real)
    (μ, A, ω) = p(t)
    dX1 = X[2]
    dX2 = μ.*(1-X[1].^2).*X[2] - X[1] + A*sin(ω*t)
    return SVector{2}(dX1, dX2)
end


"""Forcev Van der Pol Oscillator"""
function forcedvanderpol(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

forcedvanderpolSim = Process(
    process = forcedvanderpol,
    X0 = [1.0, 1.0],
    parameter_profile = (constant, constant, constant),
    parameter_profile_parameters = ((1.0,), (1.0,), (2.0,)),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.1,
    tmax = 100.0,
    alg = AutoVern9(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-9, :abs_tol => 1e-9)) # the Van der Pol equation is non-stiff for low μ, but stiff for high μ
export forcedvanderpolSim








@inline @inbounds function vanderpol(X::AbstractArray, μ::Function, t::Real)
    dX1 = X[2]
    dX2 = μ(t).*(1-X[1].^2).*X[2] - X[1]
    return SVector{2}(dX1, dX2)
end

@inline @inbounds function vanderpol_J(J, X::AbstractArray, μ::Function, t::Real)
    J .=               [  0                  1;
                    -2*μ(t)*X[1]*X[2]-1      μ(t)*(1-X[1]^2)]
end

"""Van der Pol Oscillator"""
function vanderpol(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=vanderpol_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

vanderpolSim = Process(
    process = vanderpol,
    X0 = [1.0, 1.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (5,), # (threshold, baseline, stepHeight)
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.1,
    tmax = 1000.0,
    alg = AutoVern9(Rodas5())) # the Van der Pol equation is non-stiff for low μ, but stiff for high μ
export vanderpolSim



@inline @inbounds function harmonic(X::AbstractArray, ω::Function, t::Real)
    dX2 = -ω(t)^2.0*X[1]
    dX1 = X[2]
    return SVector{2}(dX1, dX2)
end

@inline @inbounds function harmonic_J(J, X::AbstractArray, ω::Function, t::Real)
    J .= [0.0 1.0; -ω(t)^2 0.0]
end

"""Harmonic Oscillator"""
function harmonic(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=harmonic_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

harmonicSim = Process(
    process = harmonic,
    X0 = [1.0, 0.0],
    parameter_profile = unitStep,
    parameter_profile_parameters = (50.0, 2π, 4π), # (threshold, baseline, stepHeight)
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 100.0,
    alg = RK4())
export harmonicSim



@inline @inbounds function pitchfork(x::Vector, p::Function, t::Real)
    (μ, α) = p(t)
    dx = μ*x[1] + α*x[1]^3
    return [dx]
end

@inline @inbounds function pitchfork_J(J, x::Vector, p::Function, t::Real)
    (μ, α) = p(t)
    J .= [μ + 3α*x[1]^2]
end

"""Normal form for a pitchfork bifurcation"""
function pitchfork(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=pitchfork_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

pitchforkSim = Process(
    process = pitchfork,
    X0 = [0.01],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = (1.0, -1.0),# Supercritical by default
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.01,
    savedt = 0.05,
    tmax = 1000.0,
    alg = RK4(),
    solver_opts = Dict(:adaptive => false))
export pitchforkSim

@inline @inbounds function skewedHarmonic(X::AbstractArray, p::Function, t::Real)
    (ω, κ) = p(t)
    dX2 = -ω^2.0*X[1] + ω^2*X[1]*exp(-κ*X[1])
    dX1 = X[2]
    return SVector{2}(dX1, dX2)
end

@inline @inbounds function skewedHarmonic_J(J, X::AbstractArray, p::Function, t::Real)
    (ω, κ) = p(t)
    J .= [0.0 1.0; -ω^2 - κ^2*exp(-κ*x) 0.0]
end

"""Skewed 'Harmonic' Oscillator"""
function skewedHarmonic(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=skewedHarmonic_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

skewedHarmonicSim = Process(
    process = skewedHarmonic,
    X0 = [1.0, 0.0],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = ((1.0π,), (10.0,)), # (threshold, baseline, stepHeight)
    transient_t0 = -10.0,
    t0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 100.0,
    alg = RK4())
export skewedHarmonicSim
