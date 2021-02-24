# ------------------------------------------------------------------------------------------------ #
#                                   A list of default simulations                                  #
# ------------------------------------------------------------------------------------------------ #


vanderpolSim = Process(
    process = vanderpol,
    X0 = [1.0, 1.0],
    parameter_function = "unitStep",
    parameter_function_parameters = (1000.0, 0.0, 20.0), # (threshold, baseline, stepHeight)
    transient = 100.0,
    tmax = 2000.0,
    dt= 0.01,
    t0 = 0.0,
    parameter_rng = nothing,
    solver_rng = nothing,
    solver = "AutoTsit5(Rosenbrock23())", # the Van der Pol equation is non-stiff for low μ, but stiff for high μ, so we use a solver that determines which regime we are in (assuming the parameter function does not cross regimes)
    relative_tolerance = 1e-6)
export vanderpolSim


henonSim = Process(
    process = henon,
    X0 = [0.0, 0.0],
    parameter_function = ("constantParameter", "constantParameter"),
    parameter_function_parameters = ((1.4,), (0.3,)),
    transient = 100.0,
    tmax = 2000.0,
    dt= 1,
    t0 = 0.0,
    parameter_rng = nothing,
    solver_rng = nothing,
    solver = "Discrete",
    relative_tolerance = 1e-6)
export henonSim
