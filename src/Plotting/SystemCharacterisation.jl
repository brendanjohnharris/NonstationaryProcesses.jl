function labeldims(n)
    if n == 1
        ("\$t\$", "\$x\$")
    elseif n == 2
        ("\$x\$", "\$y\$")
    else
        ("\$x\$", "\$y\$", "\$z\$")
    end
end

function stationarycharacteristics(P::Process; n = length(getX0(P)), labels=labeldims(n), xlabel=labels[1], ylabel=labels[2], zlabel=(length(labels) == 3 ? labels[3] : nothing), size=(400, 350))
    p1 = marginaltrajectory(P; xlabel, ylabel, zlabel, size, mainalpha=0.35)
    display(p1)
    AMI = selfmutualinfo(P)
    p2 = plot(AMI.dims[1].val, Array((AMI isa AbstractVector ? AMI : AMI[:, end:-1:1]));  legend=true,   foreground_color_legend = nothing,
    background_color_legend = nothing, legendtitle=nothing, size,
    #xticks=[0, optimize_ticks(extrema(dims(AMI, 1).val)...; k_min = 3, k_max = 5)[1]...],
    yticks=optimize_ticks(extrema(AMI)...; k_min = 3, k_max = 6)[1],
    xlim=(0, Inf), xlabel="\$\\tau\$", ylabel="AMI")
    display(p2)
    return nothing
end
export stationarycharacteristics

function plotlyapunovresponse(P::Process, p, prange; size=(400, 350), N=1000, k=length(getX0(P)), kwargs...)
    λs = lyapunovresponse(P, p, prange, N, k)
    plot(prange, λs; size, xlabel="\$\\theta\$", ylabel="\$\\lambda\$", kwargs...)
end
export plotlyapunovresponse
