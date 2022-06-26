"""FM Signal"""
function fmWave(P::Process)
    # Only parameter is the FM signal
    T = P.transient_t0:P.savedt:P.tmax
    p = parameter_functions(P)
    Δ𝑓 = 1.0
    sol = zeros(size(T))
    pint = 0.0
    for i ∈ 2:lastindex(T)
        t = T[i]
        pint += (p(t-P.savedt) + p(t))*P.savedt/2 # Crude integration, should be fine
        sol[i] = cos(2π*(t + Δ𝑓*pint))
    end
    return sol
end

fmWaveSim = Process(
    process = fmWave,
    parameter_profile = (stepNoise),
    parameter_profile_parameters = ((0.0, 100.0), 10.0, 0.5, 0.0),
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 0.01,
    tmax = 100.0)
export fmWaveSim


"""Noisy Sine"""
function noisySine(P::Process)
    seed(P.solver_rng)
    sol = [sin(t) + parameter_function(P)(t)[1]*randn() for t in P.transient_t0:P.savedt:P.tmax]
end

noisySineSim = Process(
    process = noisySine,
    X0 = [0.0],
    parameter_profile = unitStep,
    parameter_profile_parameters = (100.0, 0.0, 0.1), # (threshold, baseline, stepHeight)
    transient_t0 = -10.0,
    dt = 0.01,
    savedt = 0.01,
    tmax = 200.0,
    alg = nothing,
    solver_opts=Dict())
export noisySineSim




"""Noisy Shifty Scaly Sine"""
function noisyShiftyScalySine(P::Process)
    # Parameters (η, C, A)
    seed(P.solver_rng)
    (η, C, A) = parameter_functions(P)
    sol = [A(t)*sin(t) + η(t)*randn() + C(t) for t in P.transient_t0:P.savedt:P.tmax]
end

function noisyShiftyShcalySine(P::Process)
    # Parameters (η, C, A)
    seed(P.solver_rng)
    (η, C, A) = parameter_functions(P)
    sol = [A(t)*(sin(t) + η(t)*randn() + C(t)) for t in P.transient_t0:P.savedt:P.tmax]
end


"""Shifty Scaly Sine"""
function shcalySine(P::Process)
    seed(P.solver_rng)
    A = parameter_function(P)
    sol = [A(t)*(sin(t) + 1.0*randn() + 1.0) for t in P.transient_t0:P.savedt:P.tmax]
end

noisyShiftyScalySineSim = Process(
    process = noisyShiftyScalySine,
    X0 = [0.0],
    parameter_profile = (ramp, constantParameter, constantParameter), # (η, C, A)
    parameter_profile_parameters = ((0.001, 2.0, 0.0), (1.0,), (1.0,)),
    transient_t0 = 0.0,
    t0 = 0.0,
    dt = 0.1,
    savedt = 0.1,
    tmax = 1000.0,
    alg = nothing,
    solver_opts=Dict())
export noisyShiftyScalySineSim

noisyShiftyShcalySineSim = Process(
    process = noisyShiftyShcalySine,
    X0 = [0.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter), # (η, C, A)
    parameter_profile_parameters = ((1.0,), (1.0,), (1.0,)),
    transient_t0 = 0.0,
    t0 = 0.0,
    dt = 0.05,
    savedt = 0.05,
    tmax = 500.0,
    alg = nothing,
    solver_opts=Dict())
export noisyShiftyShcalySineSim

shcalySineSim = Process(
    process = shcalySine,
    X0 = [0.0],
    parameter_profile = ramp, # (η, C, A)
    parameter_profile_parameters = (0.001, 1.0, 0.0),
    transient_t0 = 0.0,
    t0 = 0.0,
    dt = 0.01,
    savedt = 0.01,
    tmax = 1000.0,
    alg = nothing,
    solver_opts=Dict())
export shcalySineSim
