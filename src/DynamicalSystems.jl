using .DynamicalSystems
import .DynamicalSystems: lyapunov, lyapunovspectrum

"""
Whip up a DynamicalSystem from a Process
"""
function process2ds(P::Process)
    prob = process2problem(P)
    if getalg(P) isa FunctionMap
        d = DiscreteDynamicalSystem(prob.f, prob.u0, prob.p)
    else
        d = ContinuousDynamicalSystem(prob)
    end
end
export process2ds


"""
Calculate the largest lyapunov exponent of a process
"""
function lyapunov(P::Process)
    d = process2ds(P)
    T = gettmax(P)
    Ttr = gett0(P) - gettransient_t0(P)
    lyapunov(d, T; Ttr)
end
export lyapunov

""" Calculate the lyapunov spectrum"""
function lyapunovspectrum(P::Process, N::Integer=1000, k::Integer=length(getX0(P)))
    d = process2ds(P)
    u0 = getX0(P)
    Ttr = gett0(P) - gettransient_t0(P)
    #dt = getdt(P)
    lyapunovspectrum(d, N, k; u0, Ttr)
end
export lyapunovspectrum


"""Calculate the largest lyapunov exponent as a function of a parameter"""
function lyapunovresponse(P::Process, p, prange, N=1000, k=length(getX0(P)))
    d = process2ds(P)
    λs = Array{Float64, 2}(undef, length(prange), k)
    Threads.@threads for ip ∈ 1:length(prange)
        subP = updateparam(P, p, constantParameter, prange[ip])
        if k == 1 && !(d isa DiscreteDynamicalSystem)
            λ = lyapunov(subP)
        else
            λ = lyapunovspectrum(subP, N, k)
        end
        λs[ip, :] .= λ
    end
    return λs
end
export lyapunovresponse
