using WGLMakie
using TimeseriesTools
using DifferentialEquations
using TimeseriesSurrogates
using NonstationaryProcesses
using NonstationaryProcesses.DifferentialEquationsExt
using Foresight
using FFTW
using JSServe

transient = 20000
linewidth = 0.1

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

sx = Timeseries(lorenz)
fs = rfftfreq(size(sx, 1))
ϕ = sx[:, 1] |> rfft .|> angle

_x = sx[transient:end-transient, :]
x = _x[:, 1]

alg = :relative
α⃗ = 0:0.1:0.3
α = 0.1

f = Figure(; resolution=(1980, 800))

ss = cat(TimeseriesSurrogates.surrogate.(eachcol(sx), [PartialRandomization(α, alg)])..., dims=Var)
_s = ss[transient:end-transient, :]
s = _s[:, 1]
ax = [Axis(f[1, 1]), Axis(f[2, 1]; title="Power spectrum")]
surroplot!(ax, x, s)
axislegend(ax[2])
xlims!(ax[1], (0, 2000))
xlims!(ax[2], (0, 0.15))
ax[1].title="Time series for α = $α, alg = :$alg"

# Trajectory
ax = Axis3(f[1, 2]; title="Original trajectory")
trajectory!(ax, collect.(eachcol(_x))...; colormap=:turbo, linewidth)
shadows!(ax, collect.(eachcol(_x))...; color=(:black, 0.5), linewidth)

ax = Axis3(f[2, 2]; title="Surrogate trajectory")
trajectory!(ax, collect.(eachcol(_s))...; colormap=:turbo, linewidth)
shadows!(ax, collect.(eachcol(_s))...; color=(:black, 0.5), linewidth)

# Phasogram
# ax = Foresight.PolarAxis(f[1, 3]; title="Original phasogram")
ax = Axis(f[1, 3]; title="Original phasogram", xlabel="Frequency", ylabel="Phase")
lines!(ax, fs, ϕ)

ϕ′ = ss[:, 1] |> rfft .|> angle
# ax = Foresight.PolarAxis(f[2, 3]; title="Surrogate phasogram")
ax = Axis(f[2, 3]; title="Surrogate phasogram", xlabel="Frequency", ylabel="Phase")
lines!(ax, fs, ϕ′)

colsize!(f.layout, 1, Relative(0.5))
colgap!(f.layout, 1, Relative(0.07))
colgap!(f.layout, 2, Relative(0.0))

# save(joinpath(@__DIR__, "lorenz_surrogates.pdf"), f)



# open("index.html", "w") do io
#     println(io, """
#     <html>
#         <head>
#         </head>
#         <body>
#     """)
#     JSServe.Page(exportable=true, offline=true)
#     show(io, MIME"text/html"(), f)
#     println(io, """
#         </body>
#     </html>
#     """)
# end
