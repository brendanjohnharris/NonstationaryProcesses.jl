# ------------------------------------------------------------------------------------------------ #
#                                   A list of default simulations                                  #
# ------------------------------------------------------------------------------------------------ #


vanderpolSim = Process(
    process = vanderpol,
    X0 = [1.0, 1.0],
    parameter_function = :unitStep,
    parameter_function_parameters = (1000.0, 0.0, 20.0), # (threshold, baseline, stepHeight)
    t0 = -100.0,
    tmax = 2000.0,
    dt = 0.0001,
    savedt = 0.01,
    savet0 = 0.0,
    parameter_rng = nothing,
    solver_rng = nothing,
    solver = :(AutoTsit5(Rosenbrock23()))) # the Van der Pol equation is non-stiff for low μ, but stiff for high μ, so we use a solver that determines which regime we are in (assuming the parameter function does not cross regimes)
export vanderpolSim


henonSim = Process(
    process = henon,
    X0 = [0.0, 0.0],
    parameter_function = (:constantParameter, :constantParameter),
    parameter_function_parameters = ((1.4,), (0.3,)),
    t0 = -100.0,
    tmax = 2000.0,
    dt = 1,
    savedt = 1,
    savet0 = 0.0,
    parameter_rng = nothing,
    solver_rng = nothing,
    solver = :Discrete)
export henonSim

harmonicSim = Process(
    process = harmonic,
    X0 = [1.0, 0.0],
    parameter_function = :constantParameter,
    parameter_function_parameters = (1.0,),
    t0 = -100.0,
    tmax = 2000.0,
    dt = nothing,
    savedt = 0.01,
    savet0 = 0.0,
    parameter_rng = nothing,
    solver_rng = nothing,
    solver = :(Euler()),
    solver_opts = Dict(:adaptive=>false)) # Can't use dt with an adaptive timestep solver?
export harmonicSim
