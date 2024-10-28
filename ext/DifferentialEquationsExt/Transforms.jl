using Distributions
using Tullio
using FFTW
using StatsBase


"""
Wrap scipy's inverse short time fourier transform.

    scipy.signal.istft(Zxx, fs=1.0, window='hann', nperseg=None, noverlap=None, nfft=None, input_onesided=True, boundary=True, time_axis=- 1, freq_axis=- 2)
"""
function istft(Zxx; kwargs...)
    signal.istft(Zxx; kwargs...)
end
export istft

function stft(x; kwargs...)
    signal.stft(x; kwargs...)
end
export stft

function linearnonstationaryspectrum(g::Function, Δt, T, windows::Int=20, normalise=true; doplot=false)
    #g(x, t) is a power function
    # Run a quick stft to get the right samplings and arrays
    t⃗ = 0.0:Δt:T
    x = randn(length(t⃗))
    f, t, Z = stft(x; fs=1.0 / Δt, nperseg=floor(length(t⃗) / windows) - 1)

    R = zeros(size(Z))
    Φ = rand(Uniform(0, 2π), size(Z)) # Uniform random phase

    for ti ∈ 1:size(Z, 2)
        ff = g.(f, (t[ti],))
        if isinf(ff[1])
            ff[1] = 0.0 # So we don't diverge at 0 frequency for power laws
        end
        if normalise
            E = sum(ff .* (f[2] - f[1])) # Approximate energy
            ff = ff ./ E # So the power spectrum has unit energy
        end
        R[:, ti] = sqrt.(ff) # Fourier amplitudes
    end
    if doplot
        p = heatmap(t[1:end], log10.(f[2:end]), log10.(R[2:end, :] .^ 2), xguide="t", yguide="log₁₀(f)", colorbar_title="log₁₀(S)", yticks=-3:1, framestyle=:box, clims=(-1.5, 1.5))
        display(p)
    end
    ZZ = R .* exp.(Φ .* im)
    return istft(ZZ; fs=1.0 / Δt, nperseg=floor(length(t⃗) / windows) - 1)[end]
end
export linearnonstationaryspectrum


function selfAffine(β, Δt, T, args...)
    g(x, t) = 1 / (x^(β(t)))
    linearnonstationaryspectrum(g, Δt, T, args...)
end
export selfAffine

"""
    Generate a locally self-affine time series, which has a spectrum that varies over time according to a power-law parameterised by β(t)

Similar to Theiler1992, who uses an inverse Fourier transform for random-phase (linear) surrogates from a given power spectrum. However, this process uses the inverse short time fourier transform to construct a non-stationary signal with a time-varying power spectrum.
"""
function selfAffine(P::Process)
    # Just one parameter, β
    seed(P.solver_rng) # Maybe this doesn't propagate to python...
    β = parameter_function(P)
    Δt = getΔt(P)
    T = gettmax(P) - gettransient_t0(P)
    return selfAffine(β, Δt, T)
end

selfAffineSim = Process(
    process=selfAffine,
    parameter_profile=ramp,
    parameter_profile_parameters=(0.0, 1.0, 0.0, 2000.0),
    transient_t0=-100.0,
    t0=0.0,
    dt=0.01,
    savedt=0.1,
    tmax=2000.0)
export selfAffineSim




"""
    fouriersurrogate(x::AbstractVector, t::AbstractVector; f::Function=x->x, g::Function=x->x; kwargs...)

Generate a non-stationary surrogate time series from an input `x` after performing manipulations in (short-time) Fourier space. The function `g(A, t)` is applied to the amplitudes of the Fourier coefficients whereas `h(𝜑, t)` is applied to the phases. After modifying the short-time fourier representation of the time series, the inverse short-time Fourier transform is used to construct the surrogate time series. `x` should be regularly sampled, and the time coordinates in `t` are crucial for correctly evaluating `g` and `h`. Optional `kwargs` are passed to python's `stft` and `istft` functions.
"""
function fouriersurrogate(x::AbstractVector, t::AbstractVector; g::Function=(𝑓, 𝑡, A) -> A, h::Function=(𝑓, 𝑡, 𝜑) -> 𝜑, nperseg=1000, kwargs...)::AbstractVector
    Δt = (t[2] - t[1])
    fs = 1 / Δt
    @assert all(diff(t, dims=1) .≈ Δt) # Check regularly sampled
    𝑓, 𝑡, 𝑍 = stft(x; fs, nperseg, kwargs...)
    𝑡 .+= t[1] # Offset the start time to match input signal
    gg .= g((𝑓,), 𝑡, (abs.(𝑍[:, j]),), (angle.(𝑍[:, j]),))
    hh .= h((𝑓,), 𝑡, (abs.(𝑍[:, j]),), (angle.(𝑍[:, j]),))
    𝑍 = gg .* cis.(hh)
    𝑡, x̂ = istft(𝑍; fs, nperseg, kwargs...)
    @assert (𝑡[2] - 𝑡[1]) == (t[2] - t[1])
    x̂ = x̂[1:length(x)] # The stft is zero-padded by python to fit evenly into the windows, so the istft is longer than the input x even though it has the same sampling frequency.
end

function fouriersurrogate(x::AbstractArray, dt; kwargs...)
    # t = timeDims(x)
    x̂ = fouriersurrogate(Array(x), dt; kwargs...)
    # DimArray(x̂, (Ti(t),))
end

export fouriersurrogate


"""
Function for randomising (from a uniform distribution between 0 and 2π) an angle with a probability of 𝜂.
"""
function corruptangle(𝜑::Number, 𝜂)
    𝜑 %= 2π
    rand() < 𝜂 ? rand() * 2π : 𝜑
end
corruptangle(𝜑::Vector, 𝜂) = corruptangle.(𝜑, (𝜂,))



"""
Set phases below a threshold to 0 (synchronised)
Should prevent any issues with the desnity of randomisation interfering with power spectrum
"""
function thresholdsynchronise(𝑓::Vector, 𝜑::Number, ν)
    𝑓 < ν ? 0.0 : 𝜑
end
thresholdsynchronise(𝑓::Vector, 𝜑::Vector, ν) = thresholdsynchronise(𝑓, 𝜑, (ν,))



# """
# Take a Process and use the dark magic to produce a corrupted version, which has an extra parameter controlling the probability of the phase of each fourier coefficient being randomised. If planning to save and load this process, an instance of it must first be loaded so that the simulating function is exported
# """
# function corruptphase(P::Process, parameter_profile=constant, parameter_profile_parameters=0.0)
#     # In this case, savedt should really be a multiple of dt
#     S = P()
#     ps = getparameter_profile_parameters(P)
#     if getparameter_profile(P) isa Function || length(getparameter_profile(P)) == 1
#         S.parameter_profile = (getparameter_profile(P), parameter_profile)
#         ps = [ps]
#     else
#         S.parameter_profile = (getparameter_profile(P)..., parameter_profile)
#     end
#     S.parameter_profile_parameters = (ps..., parameter_profile_parameters)
#     pr = string(getprocess(P))
#     fname = Symbol("phaseCorrupted"*titlecase(pr))
#     pr = Symbol(pr)
#     @eval begin # Cry me a river
#         function ($fname)(S::Process)
#             D = S()
#             D.parameter_profile = getparameter_profile(S)[1:end-1]
#             if length(D.parameter_profile) == 1
#                 D.parameter_profile = D.parameter_profile[1]
#             end
#             D.parameter_profile_parameters = getparameter_profile_parameters(S)[1:end-1]
#             if length(D.parameter_profile_parameters) == 1
#                 D.parameter_profile_parameters = D.parameter_profile_parameters[1]
#             end
#             downsample = floor(D.savedt/D.dt) |> Int
#             D.savedt = D.dt # To keep accuracy for the transforms
#             D.process = $pr
#             x = timeseries(D, transient=true)
#             𝜂 = getparameter_profile(S)[end](getparameter_profile_parameters(S)[end]...)
#             ys = [Vector(x[:, i]) for i ∈ 1:size(x, 2)]
#             x = hcat([fouriersurrogate(y, times(D, transient=true); h=(𝑓, 𝑡, 𝜑)->corruptangle(𝜑, 𝜂(𝑡)), nperseg=1000*downsample) for y ∈ ys]...)
#             x = x[1:downsample:end, :]
#         end
#         export $fname
#     end
#     S.process = @eval $fname
#     return S
# end
# export corruptphase





"""
Take a Process and synchronise phases according to a frequency cutoff ν
"""
function synchronisephase(P::Process, parameter_profile=constant, parameter_profile_parameters=0.0)
    # In this case, savedt should really be a multiple of dt
    S = P()
    ps = getparameter_profile_parameters(P)
    if getparameter_profile(P) isa Function || length(getparameter_profile(P)) == 1
        S.parameter_profile = (getparameter_profile(P), parameter_profile)
        ps = [ps]
    else
        S.parameter_profile = (getparameter_profile(P)..., parameter_profile)
    end
    S.parameter_profile_parameters = (ps..., parameter_profile_parameters)
    pr = string(getprocess(P))
    fname = Symbol("phaseSynchronised" * titlecase(pr))
    pr = Symbol(pr)
    @eval begin
        function ($fname)(S::Process)
            D = S()
            D.parameter_profile = getparameter_profile(S)[1:end-1]
            if length(D.parameter_profile) == 1
                D.parameter_profile = D.parameter_profile[1]
            end
            D.parameter_profile_parameters = getparameter_profile_parameters(S)[1:end-1]
            if length(D.parameter_profile_parameters) == 1
                D.parameter_profile_parameters = D.parameter_profile_parameters[1]
            end
            downsample = floor(D.savedt / D.dt) |> Int
            D.savedt = D.dt # To keep accuracy for the transforms
            D.process = $pr
            x = timeseries(D, transient=true)
            ν = getparameter_profile(S)[end](getparameter_profile_parameters(S)[end]...)
            ys = [Vector(x[:, i]) for i ∈ 1:size(x, 2)]
            x = hcat([fouriersurrogate(y, times(D, transient=true); h=(𝑓, 𝑡, A, 𝜑) -> thresholdsynchronise(𝑓, 𝜑, ν(𝑡)), nperseg=1000 * downsample) for y ∈ ys]...)
            x = x[1:downsample:end, :]
        end
        export $fname
    end
    S.process = @eval $fname
    return S
end
export synchronisephase




"""
Now we switch gears and go for a simpler 'windowed fourier transform'; it's just the stft, but with non-overlapping, very wide windows that will match exactly the estimation windows we use. This is basically just Carl, since we won't worry about any nonsense like discontinuities at window edges.
"""
function wft(x, t; nwindows=20)
    Δt = (t[2] - t[1])
    fs = 1 / Δt
    @assert all(diff(t, dims=1) .≈ Δt) # Check regularly sampled
    remainder = length(x) % nwindows
    remainder > 1 ? (@warn "The time series does not divide well into the number of windows supplied. The remainder is $remainder") : nothing
    𝑓 = rfftfreq((length(x) - remainder) ÷ nwindows, fs)
    𝑥, 𝑡 = [reshape(i[1:end-remainder], ((length(i) - remainder) ÷ nwindows, nwindows)) for i ∈ (x, t)]
    # Have to be careful with floating points here. Start by converting all times to integers, since we know they are equally spaced
    𝑡ᵢ = Int.(round.(𝑡 ./ Δt))
    𝑡 = mean(𝑡ᵢ, dims=1)[:] .* Δt # No (fewer?) floating point errors
    𝑍 = rfft(𝑥, 1)
    return (𝑓, 𝑡, 𝑍)
end
export wft


"""
Inverse of wft transformation. Just simple windows. If the original time series did not fit well into the number of windows, give a remainder here and it will append that many samples onto the end of the time series, which all have the same value as the last point in the reconstructed time series (shouldn't be a massive deal if the remainder is only 1 or 2 and the sampling period is small).
"""
function iwft(𝑡, 𝑍; remainder=0)
    if length(𝑡) == 1
        x̂ = irfft(𝑍, 2 * size(𝑍, 1) - 1, 1) # Not usually even in this case
        x̂ = reshape(x̂, length(x̂))
        t̂ = 1:length(x̂) # Can't do a whole lot better with one window time
    else
        nwindows = size(𝑍, 2)
        @assert nwindows == length(𝑡) # 𝑡 has window centres
        Δ𝑡 = (𝑡[2] - 𝑡[1])
        @assert all(diff(𝑡, dims=1) .≈ Δ𝑡) # Equally spaced window centres
        x̂ = irfft(𝑍, 2 * size(𝑍, 1) - 2, 1)
        x̂ = reshape(x̂, length(x̂))
        Δt̂ = Δ𝑡 / ((length(x̂)) ÷ nwindows)
        N = Δ𝑡 / Δt̂
        t̂₀ = 𝑡[1] - Δt̂ * (N - 1) / 2
        append!(x̂, fill(x̂[end], remainder))
        t̂ = t̂₀:Δt̂:Δt̂*(N*nwindows-1+remainder)
    end
    return (x̂, t̂)
end
iwft(𝑓, 𝑡, 𝑍; kwargs...) = iwft(𝑡, 𝑍; kwargs...) # 𝑓 not needed, but in case you want to pass wft result directly to iwft
# iwft(𝑍::AbstractDimArray; kwargs...) = iwft(timeDims(𝑍), Array(𝑍); kwargs...)
export iwft




"""
Randomise phases of frequencies that constitute a propotion `p` of power (from the high to low frequencies)
"""
function thresholdcorrupt(𝜑, A, 𝑝) # ! SO P IS THE SQRT ROOT OF PROPORTIONAL POWER!!!
    𝜑 .= (𝜑 .+ 2π) .% 2π
    A[1] = 0 # The first amplitude is just the offset, which we don't care about (i.e. we want power of a mean-centred signal)
    psd = (A .^ 2) ./ sum(A .^ 2)
    cpsd = cumsum(psd)
    # display(plot(cpsd))
    # display(plot(A))
    idx = findfirst(cpsd .> 1 - 𝑝^2)
    isnothing(idx) ? (return 𝜑) : nothing
    phi = deepcopy(𝜑)
    phi[idx:end] .= rand(length(phi[idx:end])) .* 2π
    #display(plot([𝜑, phi]))
    return phi
end
thresholdcorrupt(𝑝) = (𝑓, 𝜑, 𝑍) -> thresholdcorrupt(𝜑, 𝑍, 𝑝)




function windowedfouriersurrogate(x::AbstractVector, Δt; g::Function=(𝑓, 𝑡, A, 𝜑) -> A, h::Function=(𝑓, 𝑡, A, 𝜑) -> 𝜑, nwindows=20)::AbstractVector
    fs = 1 / Δt
    @assert all(diff(t, dims=1) .≈ Δt) # Check regularly sampled
    𝑓, 𝑡, 𝑍 = wft(x, t; nwindows)
    #display(heatmap(log10.(abs.(𝑍[2:end-1, :])), scale=:log))
    @tullio gg[j] := g(𝑓, 𝑡[j], abs.(𝑍[:, j]), angle.(𝑍[:, j]))
    @tullio hh[j] := h(𝑓, 𝑡[j], abs.(𝑍[:, j]), angle.(𝑍[:, j]))
    gg = hcat(gg...)
    hh = hcat(hh...)
    𝑍 = gg .* cis.(hh)
    #display(plot(angle.(𝑍)))
    #@tullio 𝑍[i, j] = g(𝑓[i], 𝑡[j], abs(𝑍[i, j]), angle(𝑍[i, j]))*cis(h(𝑓[i], 𝑡[j], abs(𝑍[i, j]), angle(𝑍[i, j])))
    #display(heatmap(log10.(abs.(𝑍[2:end-1, :])), scale=:log))
    remainder = Int(length(x) % nwindows)
    x̂, 𝑡 = iwft(𝑡, 𝑍, remainder=remainder) # For most situations in which this function is a good idea, this remainder will be 1
    @assert nwindows == 1 || (𝑡[2] - 𝑡[1]) ≈ (t[2] - t[1])
    return x̂
end

function windowedfouriersurrogate(x::AbstractArray, dt; kwargs...)
    # t = timeDims(x)
    x̂ = windowedfouriersurrogate(Array(x), dt; kwargs...)
    # DimArray(x̂, (Ti(t),))
    x̂
end

export windowedfouriersurrogate


"""
Take a Process and produce a corrupted version, which has an extra parameter controlling the probability of the phase of each fourier coefficient being randomised. If planning to save and load this process, an instance of it must first be loaded so that the simulating function is exported
"""
function corruptphase(P::Process, parameter_profile=constant, parameter_profile_parameters=0.0; originalres=false, nwindows=20)
    # In this case, savedt should really be a multiple of dt
    S = P()
    ps = getparameter_profile_parameters(P)
    if getparameter_profile(P) isa Function || length(getparameter_profile(P)) == 1
        S.parameter_profile = (getparameter_profile(P), parameter_profile)
        ps = [ps]
    else
        S.parameter_profile = (getparameter_profile(P)..., parameter_profile)
    end
    S.parameter_profile_parameters = (ps..., parameter_profile_parameters)
    pr = string(getprocess(P))
    fname = Symbol("phaseCorrupted" * titlecase(pr))
    pr = Symbol(pr)
    @eval begin # Cry me a river
        function ($fname)(S::Process)
            D = S()
            D.parameter_profile = getparameter_profile(S)[1:end-1]
            if length(D.parameter_profile) == 1
                D.parameter_profile = D.parameter_profile[1]
            end
            D.parameter_profile_parameters = getparameter_profile_parameters(S)[1:end-1]
            if length(D.parameter_profile_parameters) == 1
                D.parameter_profile_parameters = D.parameter_profile_parameters[1]
            end
            if $originalres
                downsample = floor(D.savedt / D.dt) |> Int
                D.savedt = D.dt # To keep accuracy for the transforms
            else
                downsample = 1
            end
            D.process = $pr
            x = timeseries(D, transient=false)
            Y = timeseries(D, transient=true)
            𝜂 = getparameter_profile(S)[end](getparameter_profile_parameters(S)[end]...)
            ys = [Vector(x[:, i]) for i ∈ 1:size(x, 2)]
            x = hcat([windowedfouriersurrogate(y, times(D, transient=false); h=(𝑓, 𝑡, A, 𝜑) -> thresholdcorrupt(𝜑, A, 𝜂(𝑡)), nwindows=$nwindows) for y ∈ ys]...)
            idx = findfirst(D.t0 .== times(D, transient=true))
            Y[idx:end, :] = x[1:downsample:end, :]
            return Y
        end
        export $fname
    end
    S.process = @eval $fname
    return S
end
export corruptphase






"""
E.g. Lorenz attractor with corrupted phases
"""
phaseCorruptedLorenzSim = corruptphase(lorenzSim(
        X0=[0.0, -0.01, 9.0],
        parameter_profile=(constant, constant, constant),
        parameter_profile_parameters=(10.0, 28.0, 8 / 3),
        transient_t0=-100.0,
        t0=0.0,
        dt=0.001,
        savedt=0.05,
        tmax=1000.0,
        alg=AutoVern7(Rodas5()),
        solver_opts=Dict(:adaptive => true, :reltol => 1e-15)), rampInterval, (0.0, 0.25, 0.0, 1000.0))


"""
E.g. Lorenz attractor with synchronised phases
"""
phaseSynchronisedLorenzSim = synchronisephase(lorenzSim(
        X0=[0.0, -0.01, 9.0],
        parameter_profile=(constant, constant, constant),
        parameter_profile_parameters=(10.0, 28.0, 8 / 3),
        transient_t0=-100.0,
        t0=0.0,
        dt=0.001,
        savedt=0.05,
        tmax=1000.0,
        alg=AutoVern7(Rodas5()),
        solver_opts=Dict(:adaptive => true, :reltol => 1e-15)), rampInterval, (0.0, 1.0, 0.0, 1000.0))




function lowpass(𝑓, x::Vector, p)
    idxs = 𝑓 .> p
    xx = deepcopy(x)
    xx[idxs] .= 0.0
    return xx
end


"""
Pass a process through a low pass fourier filter
"""
function lowpass(P::Process, parameter_profile=constant, parameter_profile_parameters=0.0; nwindows=20)
    # In this case, savedt should really be a multiple of dt
    S = P()
    ps = getparameter_profile_parameters(P)
    if getparameter_profile(P) isa Function || length(getparameter_profile(P)) == 1
        S.parameter_profile = (getparameter_profile(P), parameter_profile)
        ps = [ps]
    else
        S.parameter_profile = (getparameter_profile(P)..., parameter_profile)
    end
    S.parameter_profile_parameters = (ps..., parameter_profile_parameters)
    pr = string(getprocess(P))
    fname = Symbol("lowpass" * titlecase(pr))
    pr = Symbol(pr)
    @eval begin
        function ($fname)(S::Process)
            D = S()
            D.parameter_profile = getparameter_profile(S)[1:end-1]
            if length(D.parameter_profile) == 1
                D.parameter_profile = D.parameter_profile[1]
            end
            D.parameter_profile_parameters = getparameter_profile_parameters(S)[1:end-1]
            if length(D.parameter_profile_parameters) == 1
                D.parameter_profile_parameters = D.parameter_profile_parameters[1]
            end
            D.process = $pr
            x = timeseries(D, transient=false)
            Y = timeseries(D, transient=true)
            𝜂 = getparameter_profile(S)[end](getparameter_profile_parameters(S)[end]...)
            ys = [Vector(x[:, i]) for i ∈ 1:size(x, 2)]
            x = hcat([windowedfouriersurrogate(y, times(D, transient=false); g=(𝑓, 𝑡, A, 𝜑) -> lowpass(𝑓, A, 𝜂(𝑡)), h=(𝑓, 𝑡, A, 𝜑) -> lowpass(𝑓, 𝜑, 𝜂(𝑡)), nwindows=$nwindows) for y ∈ ys]...)
            idx = findfirst(D.t0 .== times(D, transient=true))
            Y[idx:end, :] = x
            return Y
        end
        export $fname
    end
    S.process = @eval $fname
    return S
end
export lowpass
