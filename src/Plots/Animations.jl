
# Plot a trajectory with a trail
@recipe function f(::Type{Val{:trail}}, x, y, z; trailcolor=nothing, symmetric=false)
    if length(x) > 1
        @series begin
            if !isnothing(trailcolor)
                seriescolor := trailcolor
            end
            seriescolor --> :black
            label := nothing
            seriestype := :path
            linealpha := (symmetric ? vcat(LinRange(0, 1, length(x)/2 |> floor |> Int), LinRange(1, 0, length(x)/2 |> ceil |> Int)) : LinRange(0, 1, length(x)))
            linewidth --> 2.5
            linewidth := LinRange(0.0, plotattributes[:linewidth], length(x))
            x := x
            y := y
            z := z
        end
    end
    seriestype := :scatter
    x := x[end]
    label --> nothing
    seriescolor --> :black
    markersize --> 5
    !isnothing(y) ? y := y[end] : nothing
    !isnothing(z) ? z := y[end] : nothing
    ()
end



function movie(P::Process; vars=1:length(P.X0), downsample=1, trail=10, seriestype=:trail, kwargs...)
    X = timeseries(P, vars)
    bounds = mapslices(extrema, X)
    𝐭 = times(P)
    I = length(𝐭)
    plot(;kwargs...)
    anim = @animate for i ∈ 1:downsample:I
        print("\u1b[1GFrame $(i÷downsample)/$(I÷downsample)")
        if seriestype == :trail
            trailIdxs = (i-trail:1:i) #i-(trail*size(X, 1)÷downsample):i
            trailIdxs = trailIdxs[trailIdxs .> 0]
            Xtrail = X[trailIdxs, 1]
            Ytrail = X[trailIdxs, 2]
            plot(Xtrail, Ytrail, seriestype=:trail, markersize=20,
            color=:black,
            size=(500, 500),
            xlims=bounds[1].*1.1,
            ylims=bounds[2].*1.1,
            border=:none,
            label=:none,
            kwargs...)
        else
            bounds = widen.(mapslices(extrema, X), (0.3,))
            x = [X[1:i, v] for v in 1:lastindex(X, 2)]
            if length(bounds) < 3
                bounds = hcat(bounds, [nothing])
            end
            plot!(x..., seriestype=seriestype, xlims=bounds[1], ylims=bounds[2], zlims=bounds[3])
        end
    end every 1
    anim
end
export movie


# Bring to the boil... Can extract trail plot function to remove
function doublePendulumMovie(P::Process; downsample=1, trail=10, barcolor=:gray, ballcolor=[:black, :black], trailcolor=[:cornflowerblue, :crimson], trailwidth=2, kwargs...)
    X = timeseries(P, 1:2)
    X2 = timeseries(P, 3:4)
    bound = 1.1*max(abs.([(mapslices(extrema, X2)...)...])...)
    𝐭 = times(P)
    I = length(𝐭)
    anim = @animate for i ∈ 1:downsample:I
        print("\u1b[1GFrame $(i÷downsample)/$(I÷downsample)")

        # Bars
        plot([0, X[i, 1]], [0, X[i, 2]], seriestype=:line, linewidth=5, color=barcolor, label=:none, linealpha=0.5, xlims=(-bound, bound), ylims=(-bound, bound))
        plot!([0], [0]; seriestype=:scatter,
            markersize=20,
            markercolor=:black,
            markeralpha=0.0,
            label=:none
        )
        plot!([X[i, 1], X2[i, 1]], [X[i, 2], X2[i, 2]], seriestype=:line, linewidth=5, color=barcolor, label=:none, linealpha=0.5)

        # Ball 1
        trailIdxs = i-(trail*size(X, 1)÷downsample):i
        trailIdxs = trailIdxs[trailIdxs .> 0]
        Xtrail = X[trailIdxs, 1]
        Ytrail = X[trailIdxs, 2]
        if length(trailIdxs) > 1
            plot!(Xtrail, Ytrail, seriestype=:path,
                    color=trailcolor[1],
                    linealpha=LinRange(0, 1, length(trailIdxs)),
                    label=:none,
                    linewidth=trailwidth,
            )
        end
        # Ball 2
        Xtrail = X2[trailIdxs, 1]
        Ytrail = X2[trailIdxs, 2]
        if length(trailIdxs) > 1
            plot!(Xtrail, Ytrail, seriestype=:path,
                    color=trailcolor[2],
                    linealpha=LinRange(0, 1, length(trailIdxs)),
                    label=:none,
                    linewidth=trailwidth,
            )
        end
        plot!([X[i, 1]], [X[i, 2]]; seriestype=:scatter,
            markersize=20,
            markercolor=ballcolor[1],
            size=(1000, 1000),
            border=:none,
            label=:none,
            kwargs...
        )
        plot!([X2[i, 1]], [X2[i, 2]]; seriestype=:scatter,
            markersize=20,
            markercolor=ballcolor[2],
            size=(1000, 1000),
            border=:none,
            label=:none,
            kwargs...
        )
    end every 1
    print("\n")
    anim

end
export doublePendulumMovie
