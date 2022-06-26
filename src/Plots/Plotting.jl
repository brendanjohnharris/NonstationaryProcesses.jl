using StatsPlots
using StatsPlots.Plots
using StatsPlots.Plots.PlotMeasures
using KernelDensity

function widen(x::Union{Tuple, Vector}, s::Number=0.1)
    #μ = (x[1] + x[2])/2
    w = x[2] - x[1]
    y = (x[1] - s*w, x[2] + s*w)
end
export widen
# ------------------------------------------------------------------------------------------------ #
#                          Plot recipe for Process types (e.g. label axes)                         #
# ------------------------------------------------------------------------------------------------ #

@recipe function f(P::Process; vars=1:size(P.X0)[1], transient=false, downsample=1, colormode=nothing)
    cmode = colormode
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
        if cmode != :velocity
            seriestype --> :scatter
            markersize --> 1
            markerstrokewidth --> 0
        end
        label --> nothing
    end
    if downsample != 1
        x = tuple([xx[1:downsample:end] for xx in x]...)
    end

    if cmode == :velocity
        vx = deepcopy(timeseries(P, transient=transient))
        vx = Tuple([vx[:, i] for i in 1:size(vx)[2]])
        velocity = sqrt.(sum([(r[2:end] .- collect(r[1:end-1])).^2 for r ∈ [vx[i] for i in 1:length(vx)]], dims=1)[1]) # So colours by velocity over all vars
        seriestype := :path
        line_z := velocity
        linealpha --> 0.1
        linewidth --> 1.0
        colorbar --> nothing
    elseif cmode == :fourth
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
    elseif !isnothing(cmode) && cmode != :none
        @error "Not a supported cmode"
    end

    return x
end

include("Animations.jl")
include("Trajectories.jl")
include("Spectra.jl")
include("AMI.jl")
include("SystemCharacterisation.jl")
