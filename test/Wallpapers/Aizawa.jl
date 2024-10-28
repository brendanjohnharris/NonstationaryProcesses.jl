# # Pkg.add("NonstationaryProcessesBase")
# Pkg.add(url="https://github.com/brendanjohnharris/NonstationaryProcessesBase.jl", rev="main")
# Pkg.add(url="https://github.com/brendanjohnharris/NonstationaryProcesses.jl", rev="main")
# Pkg.add("Plots")
# Pkg.add("StatsPlots")
# Pkg.add("DifferentialEquations")
# using Plots
# using StatsPlots
using CairoMakie
using TimeseriesTools
using Foresight
set_theme!(foresight(:transparent))
using DifferentialEquations
using NonstationaryProcesses
import NonstationaryProcesses.DifferentialEquationsExt.aizawa
tmax = 15000.0
𝛿() = t -> 3.5 .+ 0.75 * sin.(2π * t ./ tmax)

S = Process(;
    process=aizawa,
    X0=[-0.78450179, -0.62887672, -0.17620268],
    parameter_profile=([constantParameter for _ in 1:3]..., 𝛿),
    parameter_profile_parameters=(0.95, 0.7, 0.6, ()),
    transient_t0=-100.0,
    t0=0.0,
    dt=0.001,
    savedt=0.001,
    tmax,
    alg=AutoVern9(Rodas5()),
    solver_opts=Dict(:adaptive => true, :reltol => 1e-11, :maxiters => 1e12))


begin
    X = TimeseriesTools.TimeSeries(S) |> eachcol .|> collect
    f = Figure(size=(5760, 5760), backgroundcolor=:transparent)
    ax = Axis3(f[1, 1]; aspect=(1, 1, 1), azimuth=π / 8, elevation=π / 5)
    hidedecorations!(ax)
    hidespines!(ax)
    trajectory!(ax, X...;
        colormode=:velocity,
        linewidth=2,
        colormap=seethrough(reverse(cgrad(:inferno)), 0.05, 0.12))
    f
    save("Aizawa.png", f)
end
