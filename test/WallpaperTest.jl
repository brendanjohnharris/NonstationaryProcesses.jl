using NonstationaryProcesses
using Plots

colorgradient = :turbo
dir = "../ProcessImages/"
mkpath(dir)

# Dict of system => variables to plot
𝒮 = Dict(simplestChaoticFlowArt()               => 2:3,
        thomasCyclicallySymmetricArt()          => 1:2,
        doubleScrollArt()                       => 1:3,
        diffusionlessLorenzArt()                => 1:3,
        piecewiseLinearHyperchaosArt()          => [1, 3],
        simplifiedLorenz4DArt()                 => [3, 2, 4],
)

for (S, vars) ∈ 𝒮
        # Square ones
        mkpath(joinpath(dir, "Square"))
        p = plot(S; vars, colormethod=:velocity, linecolor=colorgradient, background_color=:black, axis=nothing, framestyle=:none, size=(900, 900), aspect_ratio=:equal, margin=15Plots.mm, dpi=1000, linewidth=1.0);
        savefig(p, joinpath(dir, "Square", string(getprocess(S))*".png"))

        # 16:9
        mkpath(joinpath(dir, "WallpaperSquare"))
        p = plot(S; vars, colormethod=:velocity, linecolor=colorgradient, background_color=:black, axis=nothing, framestyle=:none, size=(1600, 900), aspect_ratio=:equal, margin=10Plots.mm, dpi=1000, linewidth=1.5);
        savefig(p, joinpath(dir, "WallpaperSquare", string(getprocess(S))*".png"))

        # 16:9 fill
        mkpath(joinpath(dir, "WallpaperFill"))
        p = plot(S; vars, colormethod=:velocity, linecolor=colorgradient, background_color=:black, axis=nothing, framestyle=:none, size=(1600, 900), left_margin=40Plots.mm, right_margin=40Plots.mm, dpi=1000, linewidth=1.5);
        savefig(p, joinpath(dir, "WallpaperFill", string(getprocess(S))*".png"))
end
