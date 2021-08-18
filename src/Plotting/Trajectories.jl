

# ------------------------------------------------------------------------------------------------ #
#                                     Time series distributions                                    #
# ------------------------------------------------------------------------------------------------ #
@shorthands tsdensity
@recipe function f(::Type{Val{:tsdensity}}, plt::AbstractPlot; colordensity=false)
    x, y = plotattributes[:x], plotattributes[:y]
    i = isfinite.(x) .& isfinite.(y)
    x, y = x[i], y[i]
    xlims, ylims = extrema(x), widen(extrema(y), 0.1)
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
        xguide := ""
        xticks := nothing
        yticks := nothing
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
@recipe function f(::Type{Val{:marginaltrajectory2}}, plt::AbstractPlot; marginaldownsample=1, marginallength=2000, marginalcolor=:black, colormode=:velocity, mainalpha=0.5, trajectorycolor=:black)
    cmode = colormode
    x, y = plotattributes[:x], plotattributes[:y]
    i = isfinite.(x) .& isfinite.(y)
    x, y = x[i], y[i]
    xlims, ylims = (extrema(x), extrema(y))
    top_margin --> 20Plots.mm
    right_margin --> 20Plots.mm
    legend --> :none
    framestyle --> :box

    if cmode == :velocity
        velocity = sqrt.(sum([(r[2:end] .- collect(r[1:end-1])).^2 for r ∈ [x, y]], dims=1)[1])
    elseif !isnothing(cmode) && cmode != :none
        @error "Not a supported cmode"
    end


    inset_subplots := [(1, bbox(0.0, -0.1, 1.0, 0.1)), (1, bbox(1.0, 0.0, 0.1, 1.0))]
    @series begin
        seriestype := :path
        seriescolor := marginalcolor
        ticks := nothing
        linewidth --> 2
        subplot := 2
        framestyle := :none
        x := (1:marginaldownsample:length(y))[1:marginallength]
        y := (y[1:marginaldownsample:length(y)])[1:marginallength]
    end

    @series begin
        seriestype := :path
        seriescolor := marginalcolor
        ticks := nothing
        linewidth --> 2
        subplot := 3
        framestyle := :none
        y := (length(x):-marginaldownsample:1)[1:marginallength]
        x := (x[1:marginaldownsample:length(x)])[1:marginallength]
    end

    @series begin
        seriestype := :path
        if cmode == :velocity
            line_z := velocity
        end
        linecolor := trajectorycolor
        linealpha := mainalpha
        linewidth --> 0.5
        xlims --> xlims
        ylims --> ylims
        x := x
        y := y
    end
end




@shorthands marginaltrajectory3
@recipe function f(::Type{Val{:marginaltrajectory3}}, plt::AbstractPlot; buffer=0.5, linewidth=1.0, colormode=:velocity, marginalcolor=:black, mainalpha=0.75, trajectorycolor=:black)
    cmode = colormode
    x, y, z = plotattributes[:x], plotattributes[:y], plotattributes[:z]
    i = isfinite.(x) .& isfinite.(y) .& isfinite.(z)
    x, y, z = x[i], y[i], z[i]
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
    if cmode == :velocity
        velocity = sqrt.(sum([(r[2:end] .- collect(r[1:end-1])).^2 for r ∈ [x, y, z]], dims=1)[1])
    elseif !isnothing(cmode) && cmode != :none
        @error "Not a supported cmode"
    end
    @series begin
        if isnothing(marginalcolor) && cmode == :velocity
            line_z := velocity
        end
        seriestype := :path
        linecolor := marginalcolor
        linealpha := 0.1
        linewidth := 0.5*linewidth
        x := fill(xlims[1], length(x))
        y := y
        z := z
    end

    @series begin
        if isnothing(marginalcolor) && cmode == :velocity
            line_z := velocity
        end
        seriestype := :path
        linecolor := marginalcolor
        linealpha := 0.1
        linewidth := 0.5*linewidth
        x := x
        y := fill(ylims[2], length(y))
        z := z
    end

    @series begin
        if isnothing(marginalcolor) && cmode == :velocity
            line_z := velocity
        end
        linecolor := marginalcolor
        seriestype := :path
        linealpha := 0.1
        linewidth := 0.5*linewidth
        x := x
        y := y
        z := fill(zlims[1], length(z))
    end

    @series begin
        if cmode == :velocity
            line_z := velocity
        end
        seriescolor := trajectorycolor
        seriestype := :path
        linealpha := mainalpha
        linewidth --> linewidth
        x := x
        y := y
        z := z
    end

end


function marginaltrajectory(S::Process; vars=1:length(getX0(S)), marginaldownsample=1, marginallength=2500, trajectorycolor=:turbo, marginalcolor=:black, size=(700, 600), kwargs...)
    if length(vars) == 1
        t = times(S)
        x = timeseries(S, vars)
        t = t[1:marginaldownsample:end][1:marginallength]
        x = x[1:marginaldownsample:end][1:marginallength]
        plot(t, x; xlabel="t", ylabel="x", st=:tsdensity, framestyle=nothing, xticks = optimize_ticks(extrema(t)...; k_min = 3, k_max = 5)[1], yticks = optimize_ticks(extrema(x)...; k_min = 3, k_max = 5)[1], kwargs...)
    elseif length(vars) == 2 # 2D
        marginaltrajectory2(S; vars, marginaldownsample, marginallength, trajectorycolor, marginalcolor, size,
        xticks=optimize_ticks(extrema(timeseries(S, 1))...; k_min = 3, k_max = 5)[1],
        yticks=optimize_ticks(extrema(timeseries(S, 2))...; k_min = 3, k_max = 5)[1],
        kwargs...)
    else # 3D
        marginaltrajectory3(S; vars, marginaldownsample, marginallength, mainalpha=0.25, linecolor=trajectorycolor, marginalcolor, size,
        xticks=optimize_ticks(widen(extrema(timeseries(S, 1)), 0.5)...; k_min = 3, k_max = 4)[1],
        yticks=optimize_ticks(widen(extrema(timeseries(S, 2)), 0.5)...; k_min = 3, k_max = 4)[1],
        zticks=optimize_ticks(widen(extrema(timeseries(S, 3)), 0.5)...; k_min = 3, k_max = 4)[1],
        kwargs...)
    end
end
export marginaltrajectory
