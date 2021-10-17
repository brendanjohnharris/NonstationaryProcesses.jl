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
    tmax = 8000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

𝔩𝔬𝔯𝔢𝔫𝔷 = corruptphase(lorenz, rampInterval, (0.0, 0.2, 1000.0, 8000.0), nwindows = 160)()

a = animatespectrum(𝔩𝔬𝔯𝔢𝔫𝔷, downsample=500, trail=5000, colorgradient=cgrad([:black, :crimson]), phasogram=true, nwindows=160, dpi=100)
gif(a, "./thresholdcorruptlorenzsym.mp4", fps=48)
