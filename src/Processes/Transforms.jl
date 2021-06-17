using PyCall
using Distributions

"""
Wrap scipy's inverse short time fourier transform.

    scipy.signal.istft(Zxx, fs=1.0, window='hann', nperseg=None, noverlap=None, nfft=None, input_onesided=True, boundary=True, time_axis=- 1, freq_axis=- 2)[source]¶
"""
function istft(Zxx; kwargs...)
    signal.istft(Zxx; kwargs...)
end
export istft

function stft(x; kwargs...)
    signal.stft(x; kwargs...)
end
export stft

function nonstationarysurrogate(g::Function, Δt, T, windows::Int=20, normalise=true; doplot=false)
    #g(x, t) is a power function
    # Run a 'quick' stft to get the right samplings and arrays
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
export nonstationarysurrogate


function selfAffine(β, Δt, T, args...)
    g(x, t) = 1/(x^(β(t)))
    nonstationarysurrogate(g, Δt, T, args...)
end
export selfAffine

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