using Plots
using Plots.PlotMeasures
using StatsPlots
using KernelDensity
import PyPlot

function widen(x::Union{Tuple, Vector}, s::Number=0.1)
    #μ = (x[1] + x[2])/2
    w = x[2] - x[1]
    y = (x[1] - s*w, x[2] + s*w)
end
export widen
# ------------------------------------------------------------------------------------------------ #
#                          Plot recipe for Process types (e.g. label axes)                         #
# ------------------------------------------------------------------------------------------------ #

@recipe function f(P::Process; vars=1:size(P.X0)[1], transient=false, downsample=1, colormethod=nothing)
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
        if length(vars) > 3
            vars = vars[1:3] # Can only plot in 3d
        end
        x = deepcopy(timeseries(P, vars, transient=transient))
        x = Tuple([x[:, i] for i in 1:size(x)[2]])
        if colormethod != :velocity
            seriestype --> :scatter
            markersize --> 1
            markerstrokewidth --> 0
        end
        label --> nothing
    end
    if downsample != 1
        x = tuple([xx[1:downsample:end] for xx in x]...)
    end

    if colormethod == :velocity
        velocity = sqrt.(sum([(r[2:end] .- collect(r[1:end-1])).^2 for r ∈ [x[i] for i in 1:lastindex(getX0(P))]], dims=1)[1])
        seriestype := :path
        line_z := velocity
        linealpha --> 0.1
        linewidth --> 1.0
        colorbar --> nothing
    elseif colormethod == :fourth
        if length(vars) < 4
            @error "Cannot color by the fourth variable if there are not four variables specified"
        end
        thefourth = vec(x[4])
        vars = vars[1:3]
        seriestype := :path
        line_z := thefourth
        linealpha --> 0.1
        linewidth --> 1.0
        colorbar --> nothing
    elseif !isnothing(colormethod) && colormethod != :none
        @error "Not a supported colormethod"
    end

    return x
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
            trailIdxs = i-(trail*size(X, 1)÷downsample):i
            trailIdxs = trailIdxs[trailIdxs .> 0]
            Xtrail = X[trailIdxs, 1]
            Ytrail = X[trailIdxs, 2]
            if length(trailIdxs) > 1
                plot!(Xtrail, Ytrail, seriestype=:path,
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
            # Add 3D trail?
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
@recipe function f(::Type{Val{:tsdensity}}, plt::AbstractPlot; colordensity=false)
    x, y = plotattributes[:x], plotattributes[:y]
    i = isfinite.(x) .& isfinite.(y)
    x, y = x[i], y[i]
    xlims, ylims = (extrema(x), extrema(y))
    ylims = (ylims[1] - 0.1*(ylims[1] + ylims[2])/2, ylims[2] + 0.1*(ylims[1] + ylims[2])/2)
    sz = (800, 200)
    size --> sz
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



# ------------------------------------------------------------------------------------------------ #
#                                Plot 2D trajectory with time series                               #
# ------------------------------------------------------------------------------------------------ #
@shorthands marginaltrajectory2
@recipe function f(::Type{Val{:marginaltrajectory2}}, plt::AbstractPlot;)
    x, y = plotattributes[:x], plotattributes[:y]
    i = isfinite.(x) .& isfinite.(y)
    x, y = x[i], y[i]
    xlims, ylims = (extrema(x), extrema(y))
    top_margin --> 20Plots.mm
    right_margin --> 20Plots.mm
    legend --> :none
    framestyle --> :box

    @series begin
        seriestype := :path
        seriescolor --> :black
        xlims --> xlims
        ylims --> ylims
        x := x
        y := y
    end

    inset_subplots := [(1, bbox(0.0, -0.1, 1.0, 0.1)), (1, bbox(1.0, 0.0, 0.1, 1.0))]
    @series begin
        seriestype := :path
        seriescolor --> :black
        ticks := nothing
        subplot := 2
        framestyle := :none
        x := 1:length(y)
        y := y
    end

    @series begin
        seriestype := :path
        seriescolor --> :black
        ticks := nothing
        subplot := 3
        framestyle := :none
        y := length(x):-1:1
        x := x
    end
end



# ------------------------------------------------------------------------------------------------ #
#                                Plot 3D trajectory with time series                               #
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
#                                         Decide: 2D or 3D                                         #
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
#                              Plot 3D trajectory with phase portraits                             #
# ------------------------------------------------------------------------------------------------ #

function set_pane_color(color=(0, 0, 0), ax=PyPlot.gca())
    PyPlot.svg(true)
    ax.xaxis.set_pane_color(color)
    ax.yaxis.set_pane_color(color)
    ax.zaxis.set_pane_color(color)
    f = PyPlot.gcf()
end

@shorthands marginaltrajectory3
@recipe function f(::Type{Val{:marginaltrajectory3}}, plt::AbstractPlot; buffer=0.3, linewidth=1.0, colormethod=nothing, marginalcolor=nothing, mainalpha=0.75)
    x, y, z = plotattributes[:x], plotattributes[:y], plotattributes[:z]
    i = isfinite.(x) .& isfinite.(y) .& isfinite.(z)
    x, y, z= x[i], y[i], z[i]
    if all(haskey.((plotattributes,), (:xlims, :ylims, :zlims)))
        xlims, ylims, zlims = plotattributes[:xlims], plotattributes[:ylims], plotattributes[:zlims]
    else
        xlims, ylims , zlims = widen.([extrema(x), extrema(y), extrema(z)], (buffer,))
    end
    legend --> :none
    xlims --> xlims
    ylims --> ylims
    zlims --> zlims

    # TODO Add coloring by 4th variable

    if colormethod == :velocity
        velocity = sqrt.(sum([(r[2:end] .- collect(r[1:end-1])).^2 for r ∈ [x, y, z]], dims=1)[1])
    elseif !isnothing(colormethod) && colormethod != :none
        @error "Not a supported colormethod"
    end
    if isnothing(marginalcolor) && !isnothing(colormethod) && colormethod != :none
        line_z := velocity
    end

    @series begin
        seriestype := :path
        if !isnothing(marginalcolor)
            linecolor := marginalcolor
        else
            linecolor --> :black
        end
        linealpha := 0.1
        linewidth := 0.5*linewidth
        x := fill(xlims[1], length(x))
        y := y
        z := z
    end

    @series begin
        seriestype := :path
        if !isnothing(marginalcolor)
            linecolor := marginalcolor
        else
            linecolor --> :black
        end
        linealpha := 0.1
        linewidth := 0.5*linewidth
        x := x
        y := fill(ylims[2], length(y))
        z := z
    end

    @series begin
        seriestype := :path
        if !isnothing(marginalcolor)
            linecolor := marginalcolor
        else
            linecolor --> :black
        end
        linealpha := 0.1
        linewidth := 0.5*linewidth
        x := x
        y := y
        z := fill(zlims[1], length(z))
    end

    @series begin
        seriestype := :path
        if !isnothing(colormethod) && colormethod != :none
            line_z := velocity
        end
        linealpha := mainalpha
        linewidth --> linewidth
        x := x
        y := y
        z := z
    end

end
