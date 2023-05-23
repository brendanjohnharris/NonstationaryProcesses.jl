using NonstationaryProcesses
using Plots


tom = thomasCyclicallySymmetricSim(
    X0 = [0.1, 0.0, 0.0],
    parameter_profile = constantParameter,
    parameter_profile_parameters = (0.18,),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 2000.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

𝔱𝔬𝔪 = corruptphase(tom, ramp, (0.0, 1.0, 0.0, 2000.0))()

a = animatespectrum(𝔱𝔬𝔪, downsample=100, trail=20000, nperseg=2000)
gif(a, fps=24)