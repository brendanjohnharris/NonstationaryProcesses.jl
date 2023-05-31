using CairoMakie
import TimeseriesTools.TimeSeries
using TimeseriesTools
using DifferentialEquations
using TimeseriesSurrogates
using NonstationaryProcesses
using NonstationaryProcesses.DifferentialEquationsExt
using Foresight
using FFTW

transient = 20000
linewidth = 0.05

lorenz = lorenzSim(
    X0 = [0.0, -0.01, 9.0],
    parameter_profile = (constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters = (10.0, 28.0, 8/3), # Sprott's recomendation
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.001,
    savedt = 0.025,
    tmax = 1500.0,
    alg = AutoVern7(Rodas5()),
    solver_opts = Dict(:adaptive => true, :reltol => 1e-12))

sx = Timeseries(lorenz)
fs = rfftfreq(size(sx, 1))
ϕ = sx[:, 1] |> rfft .|> angle

_x = sx[transient:end-transient, :]
x = _x[:, 1]

algs = [PartialRandomization, RelativePartialRandomization, SpectralPartialRandomization]
α⃗ = [0, 0.1, 1]

for α in α⃗
for alg in algs
@info "α = $α, alg = $alg"
f = Figure(; resolution=(800, 1200))

ss = cat(TimeseriesSurrogates.surrogate.(eachcol(sx), [alg(α)])..., dims=Var)
_s = ss[transient:end-transient, :]
s = _s[:, 1]
ax = [Axis(f[1, 1:2]), Axis(f[2, 1:2]; title="Power spectrum")]
surroplot!(ax, x, s)
axislegend(ax[2])
xlims!(ax[1], (0, 1000))
xlims!(ax[2], (0, 0.3))
ax[1].title="Time series for α = $α, alg = $alg"
ax[1].titlesize = 25

# Trajectory
ax = Axis3(f[3, 1]; title="Original trajectory")
trajectory!(ax, collect.(eachcol(_x))...; colormap=:turbo, linewidth=linewidth*2)
ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
ax.azimuth = ax.azimuth[] + 0.25
ax.elevation = ax.elevation[] + 0.25
shadows!(ax, collect.(eachcol(_x))...; color=(:slategray, 0.5), linewidth)

ax = Axis3(f[3, 2]; title="Surrogate trajectory")
ax.azimuth = ax.azimuth[] + 0.25
ax.elevation = ax.elevation[] + 0.25
trajectory!(ax, collect.(eachcol(_s))...; colormap=:turbo, linewidth=linewidth*2)
ax.xlabelvisible = ax.ylabelvisible = ax.zlabelvisible = ax.xticksvisible = ax.yticksvisible = ax.zticksvisible = ax.xticklabelsvisible = ax.yticklabelsvisible = ax.zticklabelsvisible = false
shadows!(ax, collect.(eachcol(_s))...; color=(:slategray, 0.5), linewidth)

# Phasogram
ax = Foresight.PolarAxis(f[4, 1]; title="Original phasogram", rticklabelcolor = (:black, 0.0))
# ax = Axis(f[4, 1]; title="Original phasogram", xlabel="Frequency", ylabel="Phase")
scatter!(ax, fs, ϕ, markersize=1)

ϕ′ = ss[:, 1] |> rfft .|> angle
ax = Foresight.PolarAxis(f[4, 2]; title="Surrogate phasogram", rticklabelcolor = (:black, 0.0))
# ax = Axis(f[4, 2]; title="Surrogate phasogram", xlabel="Frequency", ylabel="Phase")
scatter!(ax, fs, ϕ′, markersize=1)

rowsize!(f.layout, 1, Relative(0.2))
rowsize!(f.layout, 2, Relative(0.2))
rowsize!(f.layout, 3, Relative(0.3))
save(joinpath(@__DIR__, "lorenz_surrogates_$(α)_$alg.png"), f)

end
end
