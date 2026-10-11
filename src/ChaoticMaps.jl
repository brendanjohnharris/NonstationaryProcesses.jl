
@inline @inbounds function henon(x, p, n) # In this case, p is a tuple of functions
    dx1 = 1.0 - p(n)[1]*x[1]^2 + x[2]
    dx2 = p(n)[2]*x[1]
    return SVector{2}(dx1, dx2)
end
# Jacobian:
@inline @inbounds function henon_J(J, x, p, n)
    J .=     [-2*p(n)[1]*x[1]      1.0;
                p(n)[2]            0.0]
end

"""
Henon map
"""
function henon(P::Process) # `henon` is out-of-place, so it needs a DiscreteProblem (as for `ikeda`), not an in-place ODEFunction
    seed(P.solver_rng)
    process2solution(P)
end

henonSim = Process(
    process = henon,
    X0 = [0.0, 0.0],
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = (1.4, 0.3),
    transient_t0 = 0,
    t0 = 1000,
    dt = 1,
    savedt = 1,
    tmax = 6000,
    alg = FunctionMap()) # The only discrete solver, pretty much
export henonSim





@inline @inbounds function logistic(x, 𝑟::Function, n)
    dx = 𝑟(n)*x*(1-x[1])
    return SVector{1}(dx)
end

"""Logistic map"""
function logistic(P::Process) # Out-of-place, so a DiscreteProblem rather than an in-place ODEFunction
    seed(P.solver_rng)
    process2solution(P)
end

logisticSim = Process(
    process = logistic,
    X0 = [0.1],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (4.0,), # (threshold, baseline, stepHeight)
    transient_t0 = 0,
    t0 = 1000,
    dt = 1,
    savedt = 1,
    tmax = 6000,
    alg = FunctionMap()) # The only discrete solver, pretty much
export logisticSim




"""
Ikeda map
"""
function ikeda(x, p, n)
    X, Y = x
    𝛼, 𝛽, 𝛾, 𝜇 = p(n)
    𝜙 = 𝛽 - 𝛼/(1 + X^2 + Y^2)
    dX = 𝛾 + 𝜇*(X*cos(𝜙) - Y*sin(𝜙))
    dY = 𝜇*(X*sin(𝜙) + Y*cos(𝜙))
    return SVector{2}(dX, dY)
end

ikedaSim = Process(
    process = ikeda(P) = process2solution(P),
    X0 = [0.0, 0.0],
    parameter_profile = ntuple(x->constantParameter, 4),
    parameter_profile_parameters = (17, 0.2, 0.9, 0.85), # [6.0, 0.4, 1.0, 0.9], # Sprott
    transient_t0 = 0,
    t0 = 1000,
    dt = 1,
    savedt = 1,
    tmax = 11000,
    alg = FunctionMap())
export ikedaSim

ikedaVis = ikedaSim(tmax = 200000)
export ikedaVis
