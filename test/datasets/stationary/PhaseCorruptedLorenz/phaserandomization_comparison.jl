using TimeseriesSurrogates
using NonstationaryProcesses
using CairoMakie

if false
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

    x = TimeSeries(lorenz)

end


function progressiverandomization(x)

end
