using NonstationaryProcesses
using CairoMakie

sim = gaussianBimodalSim(
    parameter_profile=(constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters=((5.0,), (1.0,), (0.2)),  # (μ, σ, α)
    transient_t0=0.0,
    t0=0.0,
    savedt=1,
    tmax=1000.0,
    solver_rng=rand(Int))

S = simulate(sim)

p1 = density(S |> Timeseries)

S = simulate(sim(parameter_profile_parameters=((5.0,), (1.0,), (0.2))))
p2 = density(S |> Timeseries)
