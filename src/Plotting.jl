using Plots
using Plots.PlotMeasures
using StatsPlots
using KernelDensity
# ------------------------------------------------------------------------------------------------ #
#                          Plot recipe for Process types (e.g. label axes)                         #
# ------------------------------------------------------------------------------------------------ #

@recipe function f(P::Process; vars=1:size(P.X0)[1], transient=false, downsample=1)
    linecolor --> :black
    markercolor --> :black
    if length(vars) == 1
        t = times(P, transient=transient)
        x = (t, timeseries(P, vars, transient=transient))
        seriestype --> :path
        xguide --> "t"
        yguide --> "x"
        label --> nothing
        title --> String(Symbol(P.process))
    else
        x = timeseries(P, vars, transient=transient)
        x = Tuple([x[:, i] for i in 1:size(x)[2]])
        seriestype --> :scatter
        markersize --> 1
        label --> nothing
        markerstrokewidth --> 0
    end
    if downsample != 1
        x = tuple([xx[1:downsample:end] for xx in x]...)
    end
    return x
end

function movie(P::Process; vars=1:length(P.X0), downsample=1, trail=10, kwargs...)
    X = timeseries(P, vars)
    bounds = mapslices(extrema, X)
    𝐭 = times(P)
    I = length(𝐭)
    anim = @animate for i ∈ 1:downsample:I
        print("\u1b[1GFrame $(i÷downsample)/$(I÷downsample)")
        trailIdxs = i-(trail*size(X, 1)÷downsample):i
        trailIdxs = trailIdxs[trailIdxs .> 0]
        Xtrail = X[trailIdxs, 1]
        Ytrail = X[trailIdxs, 2]
        if length(trailIdxs) > 1
            plot(Xtrail, Ytrail, seriestype=:path,
                    color=:red,
                    linealpha=LinRange(0, 1, length(trailIdxs)),
                    label=:none,
                    linewidth=2,
            )
        end
        plot!([X[i, 1]], [X[i, 2]]; seriestype=:scatter,
            markersize=20,
            color=:black,
            size=(500, 500),
            xlims=bounds[1].*1.1,
            ylims=bounds[2].*1.1,
            border=:none,
            label=:none,
            kwargs...
        )
    end every 1
    anim
end



# Bring to the boil... Can extract trail plot function to remove
function doublePendulumMovie(P::Process; downsample=1, trail=10, kwargs...)
    X = timeseries(P, 1:2)
    X2 = timeseries(P, 3:4)
    bound = 1.1*max(abs.([(mapslices(extrema, X2)...)...])...)
    𝐭 = times(P)
    I = length(𝐭)
    anim = @animate for i ∈ 1:downsample:I
        print("\u1b[1GFrame $(i÷downsample)/$(I÷downsample)")

        # Bars
        plot([0, X[i, 1]], [0, X[i, 2]], seriestype=:line, linewidth=5, color=:gray, label=:none, linealpha=0.5, xlims=(-bound, bound), ylims=(-bound, bound))
        plot!([0], [0]; seriestype=:scatter,
            markersize=20,
            markercolor=:black,
            markeralpha=0.0,
            label=:none
        )
        plot!([X[i, 1], X2[i, 1]], [X[i, 2], X2[i, 2]], seriestype=:line, linewidth=5, color=:gray, label=:none, linealpha=0.5)

        # Ball 1
        trailIdxs = i-(trail*size(X, 1)÷downsample):i
        trailIdxs = trailIdxs[trailIdxs .> 0]
        Xtrail = X[trailIdxs, 1]
        Ytrail = X[trailIdxs, 2]
        if length(trailIdxs) > 1
            plot!(Xtrail, Ytrail, seriestype=:path,
                    color=:cornflowerblue,
                    linealpha=LinRange(0, 1, length(trailIdxs)),
                    label=:none,
                    linewidth=2,
            )
        end
        plot!([X[i, 1]], [X[i, 2]]; seriestype=:scatter,
            markersize=20,
            color=:black,
            size=(1000, 1000),
            border=:none,
            label=:none,
            kwargs...
        )

        # Ball 2
        Xtrail = X2[trailIdxs, 1]
        Ytrail = X2[trailIdxs, 2]
        if length(trailIdxs) > 1
            plot!(Xtrail, Ytrail, seriestype=:path,
                    color=:red,
                    linealpha=LinRange(0, 1, length(trailIdxs)),
                    label=:none,
                    linewidth=2,
            )
        end
        plot!([X2[i, 1]], [X2[i, 2]]; seriestype=:scatter,
            markersize=20,
            color=:black,
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




# ------------------------------------------------------------------------------------------------ #
#                                     Time series distributions                                    #
# ------------------------------------------------------------------------------------------------ #
@shorthands tsdensity
@recipe function f(::Type{Val{:tsdensity}}, plt::AbstractPlot; colordensity=false, densityoffset=0)
    x, y = plotattributes[:x], plotattributes[:y]
    i = isfinite.(x) .& isfinite.(y)
    x, y = x[i], y[i]
    xlims, ylims = (extrema(x), extrema(y))
    ylims = (ylims[1] - 0.1*(ylims[1] + ylims[2])/2, ylims[2] + 0.1*(ylims[1] + ylims[2])/2)
    sz = (800, 200)
    size := sz
    legend --> false
    seriescolor --> :black
    left_margin --> 5mm
    bottom_margin --> 5mm
    right_margin --> 20mm
    top_margin --> 5mm
    link := :both
    grid --> false
    #layout --> @layout [
    #    ts{0.9w,1.0h} density{0.1w,1.0h}
    #]
    @series begin
        if colordensity
            seriestype := :tscolordensity
        else
            seriestype := :path
        end
        framestyle --> :on
        xguide --> "t"
        xlims --> xlims
        ylims --> ylims
        x := x
        y := y
    end


    idxs = LinRange(ylims..., 1000)
    p = pdf(kde(y), idxs)
    @series begin
        seriestype := :shape
        ticks := nothing
        inset_subplots := (1, bbox(1.0, 0.0, 0.1, 1.0))
        subplot := 2
        background_color_inside := nothing
        foreground_color_border := nothing
        fillcolor --> plotattributes[:seriescolor]
        linecolor --> plotattributes[:seriescolor]
        framestyle --> :none
        yguide := ""
        ylims --> ylims
        xlims --> (0, max(p...))
        fillalpha := 0.25
        x := [0.0, p..., 0.0]
        y := [ylims[1], idxs..., ylims[2]]
    end
end

tsdensity(P::Process; vars=1, kwargs...) = tsdensity(timeseries(P, vars); kwargs...)
export tsdensity

@shorthands tscolordensity
@recipe function f(::Type{Val{:tscolordensity}}, plt::AbstractPlot)
    x, y = plotattributes[:x], plotattributes[:y]
    i = isfinite.(x) .& isfinite.(y)
    x, y = x[i], y[i]
    xlims, ylims = (extrema(x), extrema(y))
    ylims = (ylims[1] - 0.1*(ylims[1] + ylims[2])/2, ylims[2] + 0.1*(ylims[1] + ylims[2])/2)
    size := (800, 200)
    legend --> false
    seriescolor --> :black
    link := :both
    bottom_margin --> 5mm
    grid --> false


    #idxs = LinRange(ylims..., 1000)
    p = pdf(kde(y), y)

    @series begin
        seriestype := :path
        right_margin --> 0mm
        top_margin --> 0mm
        colorbar --> false
        framestyle --> :on
        xguide --> "t"
        subplot := 1
        xlims --> xlims
        ylims --> ylims
        line_z --> p
        x := x
        y := y
    end

end
tscolordensity(P::Process; vars=1, kwargs...) = tscolordensity(timeseries(P, vars); kwargs...)
export tscolordensity