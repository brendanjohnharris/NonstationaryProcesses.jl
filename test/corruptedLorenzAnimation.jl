using NonstationaryProcesses
using Plots

lorenz = lorenzSim(
    X0 = [0.0, -0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3), # Sprott's recomendation
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 1000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

𝔩𝔬𝔯𝔢𝔫𝔷 = corruptphase(lorenz, rampInterval, (0.0, 0.25, 500.0, 1000.0))()

a = animatespectrum(𝔩𝔬𝔯𝔢𝔫𝔷, downsample=100, trail=5000, colorgradient=cgrad([:black, :crimson]), phasogram=true, nperseg=1000, dpi=100)
gif(a, "./thresholdcorruptlorenzsym.gif", fps=24)