# ------------------------------------------------------------------------------------------------ #
#                                   A list of default simulations                                  #
# ------------------------------------------------------------------------------------------------ #


vanderpolSim = Process(
    process = vanderpol,
    X0 = [1.0, 1.0],
    parameter_profile = unitStep,
    parameter_profile_parameters = (1000.0, 0.0, 20.0), # (threshold, baseline, stepHeight)
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.1,
    tmax = 2000.0,
    alg = RK4()) # the Van der Pol equation is non-stiff for low μ, but stiff for high μ
export vanderpolSim

henonSim = Process(
    process = henon,
    X0 = [0.0, 0.0],
    parameter_profile = (unitStep, constantParameter),
    parameter_profile_parameters = ((1000, 1.4, -0.4), (0.3,)), # (threshold, baseline, stepHeight)
    transient_t0 = 0,
    t0 = 0,
    dt = 1,
    savedt = 1,
    tmax = 2000,
    alg = FunctionMap()) # The only discrete solver, pretty much
export henonSim

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


noisyShiftyScalySineSim = Process(
    process = noisyShiftyScalySine,
    X0 = [0.0],
    parameter_profile = (constantParameter, ramp, ramp), # (η, C, A)
    parameter_profile_parameters = ((0.0), (0.001, 1.0, 0.0), (0.001, 1.0, 0.0)),
    transient_t0 = 0.0,
    t0 = 0.0,
    dt = 0.01,
    savedt = 0.01,
    tmax = 1000.0,
    alg = nothing,
    solver_opts=Dict())
export noisyShiftyScalySineSim


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


# skewedHarmonicSim = Process(
#     process = skewedHarmonic,
#     X0 = [1.0, 0.0],
#     parameter_profile = (constantParameter, constantParameter),
#     parameter_profile_parameters = ((1.0π,), (10.0,)), # (threshold, baseline, stepHeight)
#     transient_t0 = -10.0,
#     t0 = 0.0,
#     dt = 0.0001,
#     savedt = 0.001,
#     tmax = 100.0,
#     alg = RK4())
# export skewedHarmonicSim



skewedGaussianQuadraticSim = Process(
    process = skewedGaussianQuadratic,
    X0 = [0.0],
    parameter_profile = (constantParameter, constantParameter, ramp),
    parameter_profile_parameters = ((-1.0,), (2.0,), (0.01, 0.1, 0.0)),
    transient_t0 = -1.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.001,
    tmax = 1000.0,
    alg = EM())
export skewedGaussianQuadraticSim


gaussianBimodalSim = Process(
    process = gaussianBimodal,
    parameter_profile = (ramp, constantParameter, constantParameter),
    parameter_profile_parameters = ((0.001, 0.1, 0.0), (1.0,), (0.1)),
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 1,
    tmax = 10000.0)
export gaussianBimodalSim


shiftyNoiseSim = Process(
    process = shiftyNoise,
    parameter_profile = (constantParameter, stepNoise),
    parameter_profile_parameters = [(1.0,), ((0.0, 1000.0), 100.0, 2.0, 0.0)],
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 1,
    tmax = 1000.0)
export shiftyNoiseSim


fmWaveSim = Process(
    process = fmWave,
    parameter_profile = (stepNoise),
    parameter_profile_parameters = ((0.0, 100.0), 10.0, 0.5, 0.0),
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 0.01,
    tmax = 100.0)
export fmWaveSim
