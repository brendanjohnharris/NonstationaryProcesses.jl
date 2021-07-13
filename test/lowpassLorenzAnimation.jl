using NonstationaryProcesses
using Plots

savedt = 0.01
S = lorenzSim(
    X0 = [0.0, -0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3), # Sprott's recomendation
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = savedt,
    tmax = 1000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

P = NonstationaryProcesses.lowpass(S, rampInterval, (1/(2*savedt), 0.0, 200.0, 1000.0), nwindows = 100)()

a = animatespectrum(P, downsample=50, trail=5000, colorgradient=cgrad([:black, :crimson]), phasogram=true, nwindows=100, dpi=100)
gif(a, "./lowpasslorenzsym.mp4", fps=48)