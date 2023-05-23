using NonstationaryProcesses
using Plots

function wallPlot(S, vars)
        # Square ones
        if vars[2] isa String
                vars, name = vars
        else
                name = getprocess(S)
        end
        mkpath(joinpath(dir, "Square"))
        p = plot(S; vars, colormode=:velocity, linecolor=colorgradient, background_color=:black, axis=nothing, framestyle=:none, size=(900, 900), aspect_ratio=:equal, margin=15Plots.mm, dpi=1000, linewidth=1.0);
        savefig(p, joinpath(dir, "Square", string(name)*".png"))

        # 16:9
        mkpath(joinpath(dir, "WallpaperSquare"))
        p = plot(S; vars, colormode=:velocity, linecolor=colorgradient, background_color=:black, axis=nothing, framestyle=:none, size=(1600, 900), aspect_ratio=:equal, margin=10Plots.mm, dpi=1000, linewidth=1.5);
        savefig(p, joinpath(dir, "WallpaperSquare", string(name)*".png"))

        # 16:9 fill
        mkpath(joinpath(dir, "WallpaperFill"))
        p = plot(S; vars, colormode=:velocity, linecolor=colorgradient, background_color=:black, axis=nothing, framestyle=:none, size=(1600, 900), left_margin=40Plots.mm, right_margin=40Plots.mm, dpi=1000, linewidth=1.5);
        savefig(p, joinpath(dir, "WallpaperFill", string(name)*".png"))
end


colorgradient = :turbo
dir = "../ProcessImages/"
mkpath(dir)

# Dict of system => variables to plot
𝒮 = Dict(
        simplestChaoticFlowVis()                => 2:3,
        thomasCyclicallySymmetricVis()          => 1:2,
        doubleScrollVis()                       => 1:3,
        diffusionlessLorenzVis()                => 1:3,
        piecewiseLinearHyperchaosVis()          => [1, 3],
        simplifiedLorenz4DVis()                 => [3, 2, 4],
        chensSystemVis()                        => [1, 2, 3],
        chensSystemVis(parameters=(46.0, 11.0, 29.0), tmax=500.0)
                                                => ([1, 2, 3], "chensSystem46"),
        cartesianDoublePendulumSim(profiles=(constant, constant, constant, constant),
                                   parameters=(0.5, 1.0, 1.0, 1.0),
                                   X0=[π/4, π, 0.0, 0.0],
                                   tmax=2000.0,
                                   savedt= 0.01)
                                                => [3, 4],
        lorenzVis()                             => 1:3,
)



for (S, vars) ∈ 𝒮
        wallPlot(S, vars)
end
