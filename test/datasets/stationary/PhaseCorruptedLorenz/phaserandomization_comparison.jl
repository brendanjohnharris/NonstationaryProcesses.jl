using CairoMakie
using TimeseriesTools
using DifferentialEquations
using TimeseriesSurrogates
using NonstationaryProcesses
using NonstationaryProcesses.DifferentialEquationsExt


lorenz = lorenzSim(
    X0 = [0.0, -0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3), # Sprott's recomendation
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.01,
    tmax = 800.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

x = Timeseries(lorenz, 1)

algs = [:absolute, :relative]
α⃗ = 0:0.1:0.3
f = Figure(; resolution=(1980, 720))
panel = 'a' - 1
for (j, α) ∈ enumerate(α⃗)
    for (i, alg) ∈ enumerate(algs)
        global panel += 1
        @info (i, j)
        s = TimeseriesSurrogates.surrogate(x, PartialRandomization(α, alg))

        ax = [Axis(f[i, j][k, 1]) for k ∈ 1:2]
        surroplot!(ax, x[2500:end-2500], s[2500:end-2500])
        xlims!(ax[1], (0, 2000))
        xlims!(ax[2], (0, 0.15))
        ax[1].title="$(panel)): α = $α, alg = $alg"
    end
end
colgap!(f.layout, 100)
rowgap!(f.layout, 50)

save(joinpath(@__DIR__, "lorenz_surrogates.pdf"), f)
