using NonstationaryProcesses
using DynamicalSystems
using Distributions
using StatsPlots

sim = lorenzSim(
    X0 = [0.0, -0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3), # Sprott's recomendation
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.05,
    tmax = 50.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-15))

# High dim.
N = 2000
𝐃 = [Uniform(9.0, 14.0), Uniform(26.0, 56.0), Uniform(1.0, 2.8)];
ps = [rand.(𝐃) for i in 1:N];
𝒮 = [sim(ps = ps[i]) for i ∈ 1:N];

𝝀 = lyapunov.(𝒮)
density(𝝀, color=cornflowerblue, linewidth=2.5, xguide="λ", yguide="Density", label=nothing, dpi=600, bg=nothing, bginside=nothing)
