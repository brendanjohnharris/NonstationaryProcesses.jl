using PyCall
using Distributions
using Tullio


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
    f, t, Z = stft(x; fs=1.0/Δt, nperseg=floor(length(t⃗)/windows)-1)

    R = zeros(size(Z))
    Φ = rand(Uniform(0, 2π), size(Z)) # Uniform random phase

    for ti ∈ 1:size(Z, 2)
        ff = g.(f, (t[ti],))
        if isinf(ff[1])
            ff[1] = 0.0 # So we don't diverge at 0 frequency for power laws
        end
        if normalise
            E = sum(ff.*(f[2]-f[1])) # Approximate energy
            ff = ff./E # So the power spectrum has unit energy
        end
        R[:, ti] = sqrt.(ff) # Fourier amplitudes
    end
    if doplot
        p = heatmap(t[1:end], log10.(f[2:end]), log10.(R[2:end, :].^2), xguide="t", yguide="log₁₀(f)", colorbar_title="log₁₀(S)", yticks=-3:1, framestyle=:box, clims=(-1.5, 1.5))
        display(p)
    end
    ZZ = R.*exp.(Φ.*im)
    return istft(ZZ; fs=1.0/Δt, nperseg=floor(length(t⃗)/windows)-1)[end]
end
export linearnonstationaryspectrum


function selfAffine(β, Δt, T, args...)
    g(x, t) = 1/(x^(β(t)))
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
    process = selfAffine,
    parameter_profile = ramp,
    parameter_profile_parameters = (0.0, 1.0, 0.0, 2000.0),
    transient_t0 = -100.0,
    t0 = 0.0,
    dt = 0.01,
    savedt = 0.1,
    tmax = 2000.0)
export selfAffineSim




"""
    fouriersurrogate(x::AbstractVector, t::AbstractVector; f::Function=x->x, g::Function=x->x; kwargs...)

Generate a non-stationary surrogate time series from an input `x` after performing manipulations in (short-time) Fourier space. The function `g(A, t)` is applied to the amplitudes of the Fourier coefficients whereas `h(𝜑, t)` is applied to the phases. After modifying the short-time fourier representation of the time series, the inverse short-time Fourier transform is used to construct the surrogate time series. `x` should be regularly sampled, and the time coordinates in `t` are crucial for correctly evaluating `g` and `h`. Optional `kwargs` are passed to python's `stft` and `istft` functions.
"""
function fouriersurrogate(x::AbstractVector, t; g::Function=(𝑓, 𝑡, A)->A, h::Function=(𝑓, 𝑡, 𝜑)->𝜑, nperseg=1000, kwargs...)
    Δt = (t[2] - t[1])
    fs = 1/Δt
    @assert all((t[2:end] .- t[1:end-1]) .≈ Δt) # Check regularly sampled
    𝑓, 𝑡, 𝑍 = stft(x; fs=fs, nperseg=nperseg, kwargs...)
    𝑡 .+= t[1] # Offset the start time to match input signal
    @tullio 𝑍[i, j] = g(𝑓[i], 𝑡[j], abs(𝑍[i, j]))*cis(h(𝑓[i], 𝑡[j], angle(𝑍[i, j])))
    𝑡, x̂ = istft(𝑍; fs=fs, nperseg=nperseg, kwargs...)
    @assert (𝑡[2] - 𝑡[1]) == (t[2] - t[1])
    x̂ = x̂[1:length(x)] # The stft is zero-padded by python to fit evenly into the windows, so the istft is longer than the input x even though it has the same sampling frequency.
end

function fouriersurrogate(x::DimArray; kwargs...)
    t = timeDims(x)
    x̂ = fouriersurrogate(Array(x), t; kwargs...)
    DimArray(x̂, (Ti(t),))
end

export fouriersurrogate


"""
Function for randomising (from a uniform distribution between 0 and 2π) an angle with a probability of 𝜂.
"""
function corruptangle(𝜑, 𝜂)
    𝜑 %= 2π
    rand() < 𝜂 ? rand()*2π : 𝜑
end



"""
Take a Process and use the dark magic to produce a corrupted version, which has an extra parameter controlling the probability of the phase of each fourier coefficient being randomised.
"""
function corruptphase(P::Process, parameter_profile=constant, parameter_profile_parameters=0.0)
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
    fname = Symbol("phaseCorrupted"*titlecase(pr))
    pr = Symbol(pr)
    @eval begin # 😢
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
            x = timeseries(D, transient=true)
            𝜂 = getparameter_profile(S)[end](getparameter_profile_parameters(S)[end]...)
            ys = [Vector(x[:, i]) for i ∈ 1:size(x, 2)]
            x = hcat([fouriersurrogate(y, times(S, transient=true); h=(𝑓, 𝑡, 𝜑)->corruptangle(𝜑, 𝜂(𝑡))) for y ∈ ys]...)
        end
    end
    S.process = @eval $fname
    return S
end
export corruptphase




