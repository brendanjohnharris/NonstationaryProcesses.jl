using Plots
using StatsPlots
using DifferentialEquations
using NonstationaryProcesses
import NonstationaryProcesses.DifferentialEquationsExt as DE

tmax = 1000.0
𝑘() = t -> 0.55 .+ (9 .* 0.55) * sin.(4π * t ./ tmax)
S = Process(;
    process=DE.dequanLi,
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


p = plot(S;
    vars=1:3,
    colormode=:velocity,
    linecolor=:turbo,
    background_color=:black,
    axis=nothing,
    framestyle=:none,
    size=(900, 900),
    aspect_ratio=:equal,
    margin=15Plots.mm,
    dpi=1000, linewidth=1.0,
    N=1000000,
    linealpha=0.01)

savefig(p, "DequanLi.png")
