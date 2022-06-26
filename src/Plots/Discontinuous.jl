using Plots

function Plots.plot(xr, D::Discontinuous; kwargs...)
    Plots.plot(xr, D(xr); kwargs...)
end
function Plots.plot(D::Discontinuous; xguide="t", yguide="p(t)", kwargs...) # Actually this can be a recipe?
    xr = extrema(collect(D.d))
    scale = abs(-(xr...))
    if scale == 0.0
        xr = collect(D.d)[1]-1:0.001:collect(D.d)[1]+1
    else
        xr = (xr[1] - 0.2*scale):(scale/1000):(xr[2] + 0.2*scale)
    end
    Plots.plot(xr, D; xguide, yguide, kwargs...)
end
