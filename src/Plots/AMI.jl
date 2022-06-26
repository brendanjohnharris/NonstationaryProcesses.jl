using DelayEmbeddings
using DimensionalData
import DelayEmbeddings.selfmutualinfo
using Plots


function selfmutualinfo(P::Process, lags; vars=1:length(getX0(P)), kwargs...)
    dt = getsavedt(P)
    X = timeseries(P, vars)
    AMI = selfmutualinfo.(eachcol(X), (lags,); kwargs...)
    X isa AbstractVector ? DimArray(AMI[1], Dim{:τ}(lags.*dt), :AMI) : DimArray(hcat(AMI...), (Dim{:τ}(lags.*dt), dims(X, 2)), :AMI)
end
# ! The DelayEmbeddings.mincrossing isn't so greate, because any small noise will be counted as a local minimum if the sampling period is small.
function selfmutualinfo(P::Process; τmultiple=10, kwargs...) # This one automatically calculates some sensible lags and returns an AMI function over a nice range of lags
    ts = 1:length(times(P))÷10
    AMI = selfmutualinfo(P, ts; kwargs...)
    τ = round(Int, mean(DelayEmbeddings.mincrossing.(eachcol(AMI), (1:size(AMI, 1),))))
    return AMI isa AbstractVector ? AMI[1:τ*τmultiple] : AMI[1:τ*τmultiple, :]
end
export selfmutualinfo
