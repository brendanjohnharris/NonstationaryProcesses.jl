

"""
Plot and animate a system evolving alongside its power spectrum
"""
function animatespectrum(S::Process; downsample=100, trail=1000, nwindows=20, phasogram=false, colorgradient=cgrad([:black, :crimson]), dpi=200, dospectrum=true, linewidth=1.5)
    # Phasogram really requires a high sampling period
    t = times(S)
    p = parameterseries(S)[end, :]
    X = timeseries(S)
    remainder = length(t)%nwindows
    X1X, t1t = [reshape(i[1:end-remainder], ((length(i)-remainder)÷nwindows, nwindows)) for i ∈ (vec(X[:, 1]), t)]
    bounds = mapslices(extrema, X, dims=1)
    x = vec(X[:, 1])
    y = vec(X[:, 2])
    if size(X, 2) > 3
        z = vec(X[:, 3])
    else
        z = fill(nothing, size(x))
    end
    xprev = zeros(size(x))
    anim = @animate for i = 1:downsample:length(x)
            print("\u1b[1GFrame $(i÷downsample)/$(length(x)÷downsample)")
            trailidxs = i-trail:i
            trailidxs = trailidxs[trailidxs .> 0]
            if isempty(trailidxs)
                trailidxs = [1]
            end
            widx = findfirst(mean(t1t, dims=1)[:].>=t[i])
            if isnothing(widx)
                widx = size(t1t, 2)
            end
            if dospectrum
                plot(t1t[:, widx], X1X[:, widx], seriestype=:spectrogram, teval=1, layout=(phasogram ? (@layout [a{0.8h}; b{0.65w} c{0.35w}]) : (@layout [a{0.7h}; b])), subplot=2, size=(1000, 1000), framestyle=:box, left_margin=10Plots.mm, ylims=(10^-5, 10^10), nwindows=1, axis=nothing, grid=true, colorbar=false, color=:black, gridlinewidth=4.0, dpi=dpi)
            end
            if phasogram
                # Gonna be very tricky here. Only take the top half of frequencies, and perform all sorts of rotations to get a nice visualisation
                alignshift = 0.0# -mean(x[Int(ceil(length(x)/2)):end]) # Becuase only relative phases really matter, rotate so that the phases betwene frames are mostly aligned
                alignreflect = 1.0 # sum(x .- xprev) > sum(x.+alignshift) ? -1.0 : 1.0 # Should we also reflect?
                plot!(t1t[:, widx], X1X[:, widx], seriestype=:spectrogram, teval=1, subplot=3, nwindows=1, phasogram=true, colorbar=nothing, topfreqs = nothing, axis=nothing, gridlinewidth=4.0, margin = 2Plots.mm, foreground_text_color=:white)
            end
            if size(X, 2) == 2
                plot!(x[trailidxs], y[trailidxs], z[trailidxs]; seriestype=:trail, subplot=1, axis=nothing, framestyle=:none, xlims=bounds[1].+[-0.1, 0.1], ylims=bounds[2].+[-0.1, 0.1], zlims=bounds[3].+[-0.1, 0.1], markersize=0.0, linez=p[trailidxs], clims=extrema(p), colorbar=nothing, color=colorgradient, linewidth)
            else
                plot!(x[trailidxs], y[trailidxs]; seriestype=:trail, subplot=1, axis=nothing, framestyle=:none, xlims=bounds[1].+[-0.1, 0.1], ylims=bounds[2].+[-0.1, 0.1], markersize=0.0, linez=p[trailidxs], clims=extrema(p), colorbar=nothing, color=colorgradient, symmetric=true, linewidth)
            end
            xprev = copy(x)
    end
    anim
end
export animatespectrum


# Plot the spectrogram of a time series as a heatmap. Optionally give a time at which to evaluate the spectrogram, in which case the plot will be a path.
@shorthands spectrogram
@recipe function f(::Type{Val{:spectrogram}}, x, y, z; teval=nothing, nwindows=20, nperseg=length(x)÷nwindows, phasogram=false, topfreqs=nothing)
    t = x
    Δt = (t[2] - t[1])
    fs = 1/Δt
    #@assert all((t[2:end] .- t[1:end-1]) .≈ Δt) # Check regularly sampled
    #𝑓, 𝑡, 𝑍 = stft(y; fs=fs, nperseg=nperseg)
    𝑓, 𝑡, 𝑍 = wft(y, t; nwindows)
    𝑡 .+= t[1]
    if phasogram # Show phases instead of spectrum
        S = angle.(𝑍)
        if !isnothing(topfreqs)
            idxs = Int(ceil(length(𝑓)/topfreqs)):length(𝑓)
            𝑓 = 𝑓[idxs]
            S = S[idxs, :]
        end
    else
        S = abs.(𝑍).^2
    end
    if isnothing(teval)
        seriestype := :heatmap
        if phasogram
            scale := :linear
            colorbar_title := "ϕ"
            z := Surface((S[2:end, 1:end-1]))
        else
            colorbar_title := "log₁₀(S)"
            z := Surface(log10.(S[2:end, 1:end-1]))
        end
        yscale --> :log
        x := 𝑡[1:end-1]
        y := 𝑓[2:end]
    else
        seriestype := (phasogram ? :phasogram : :path)
        _, tidx = findmin(abs.(𝑡.-teval))
        if phasogram
            title --> ""
            dolog := false
        else
            scale --> :log
            xguide --> "f"
            yguide --> "S"
            linewidth --> 2.5
        end
        legend --> nothing
        seriescolor --> :black
        x := 𝑓[2:end]
        y := S[2:end, tidx]
    end
    ()
end



@recipe function f(::Type{Val{:phasogram}}, x, y, z; dolog=false)
    projection := :polar
    seriestype := :scatter
    x := dolog ? y[x .> 0] : y
    y := dolog ? log10.(x[x .> 0]) .- min(log10.(x[x .> 0])...) : x
    markersize --> 2
    markercolor --> :black
    title --> (dolog ? "log₁₀(f)" : "f") #(dolog ? "log₁₀(f) + min[log₁₀(f)]" : "f")
    label --> nothing
    ()
end
