# ------------------------------------------------------------------------------------------------ #
#                                   A list of default simulations                                  #
# ------------------------------------------------------------------------------------------------ #


vanderpolSim = Process(
    process = vanderpol,
    X0 = [1.0, 1.0],
    parameter_profile = unitStep,
    parameter_profile_parameters = (1000.0, 0.0, 20.0), # (threshold, baseline, stepHeight)
    t0 = -100.0,
    savet0 = 0.0,
    dt = 0.001,
    savedt = 0.1,
    tmax = 2000.0,
    alg = AutoTsit5(Rosenbrock23())) # the Van der Pol equation is non-stiff for low μ, but stiff for high μ
export vanderpolSim

henonSim = Process(
    process = henon,
    X0 = [0.0, 0.0],
    parameter_profile = (unitStep, constantParameter),
    parameter_profile_parameters = ((1000, 1.4, -0.4), (0.3,)), # (threshold, baseline, stepHeight)
    t0 = 0,
    savet0 = 0,
    dt = 1,
    savedt = 1,
    tmax = 2000,
    alg = FunctionMap())
export henonSim

harmonicSim = Process(
    process = harmonic,
    X0 = [1.0, 0.0],
    parameter_profile = unitStep,
    parameter_profile_parameters = (100.0, 1.5, 1.5), # (threshold, baseline, stepHeight)
    t0 = -10.0,
    savet0 = 0.0,
    dt = 0.0001,
    savedt = 0.001,
    tmax = 200.0,
    alg = AutoTsit5(Rosenbrock23()))
export harmonicSim

