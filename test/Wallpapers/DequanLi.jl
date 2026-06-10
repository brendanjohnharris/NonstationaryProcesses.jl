# # Pkg.add("NonstationaryProcessesBase")
# Pkg.add(url="https://github.com/brendanjohnharris/NonstationaryProcessesBase.jl", rev="main")
# Pkg.add(url="https://github.com/brendanjohnharris/NonstationaryProcesses.jl", rev="main")
# Pkg.add("Plots")
# Pkg.add("StatsPlots")
# Pkg.add("DifferentialEquations")
using CairoMakie
using TimeseriesBase
using Foresight
set_theme!(foresight(:transparent))
using DifferentialEquations
using NonstationaryProcesses
import NonstationaryProcesses.dequanLi
tmax = 2000.0
𝑘() = t -> 0.55 .+ (9 .* 0.55) * sin.(4π * t ./ tmax)

S = Process(;
    process=dequanLi,
    X0=[-11.395722, -104.12614, 202.46173],
    parameter_profile=Tuple([[constantParameter for _ in 1:5]..., 𝑘]),
    parameter_profile_parameters=(40, 1.833, 0.16, 0.65, 20, ()), # (𝑎, 𝑐, 𝑑, 𝑒, 𝑓, 𝑘)
    transient_t0=-100.0,
    t0=0.0,
    dt=0.001,
    savedt=0.001,
    tmax,
    alg=AutoVern9(Rodas5()),
    solver_opts=Dict(:adaptive => true, :reltol => 1e-11, :maxiters => 1e12))

# begin
#     p = plot(S;
#         vars=1:3,
#         colormode=:velocity,
#         linecolor=seethrough(:turbo, 0.05, 0.1),
#         background_color=:white,
#         axis=nothing,
#         framestyle=:none,
#         size=(900, 900),
#         aspect_ratio=:equal,
#         margin=0Plots.mm,
#         dpi=1000,
#         linewidth=0.1,
#         N=6000000,
#         linealpha=0.05)

#     savefig(p, "DequanLi.png")
# end


begin
    X = TimeseriesBase.Timeseries(S) |> eachcol .|> collect
    f = Figure(resolution=(7680, 7680), backgroundcolor=:transparent)
    ax = Axis3(f[1, 1]; aspect=(1, 1, 1), azimuth=-π / 10, elevation=π / 5)
    hidedecorations!(ax)
    hidespines!(ax)
    trajectory!(ax, X...;
        colormode=:velocity,
        linewidth=3,
        colormap=seethrough(:viridis, 0.04, 0.1))
    f
    save("DequanLi.png", f)
end
