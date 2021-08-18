using .DynamicalSystems


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
function DynamicalSystems.lyapunov(P::Process)
    d = process2ds(P)
    T = gettmax(P)
    lyapunov(d, T)
end
export lyapunov

""" Calculate the lyapunov spectrum"""
function DynamicalSystems.lyapunovspectrum(P::Process, N=10000)
    d = process2ds(P)
    u0 = getX0(P)
    Ttr = gett0(P) - gettransient_t0(P)
    dt = getdt(P)
    lyapunovspectrum(d, N; u0, Ttr, dt)
end
export lyapunovspectrum


"""Calculate the largest lyapunov exponent as a function of a parameter"""
function lyapunovresponse(P::Process, p, prange)
    d = process2ds(P)
    λs = Array{Float64, 2}(undef, length(prange), dimension(d))
    Threads.@threads for ip ∈ 1:length(prange)
        subP = updateparam(P, p, constantParameter, prange[ip])
        λs[ip, :] = lyapunovspectrum(subP)
    end
    return λs
end
export lyapunovresponse
