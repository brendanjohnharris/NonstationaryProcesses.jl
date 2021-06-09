using .DynamicalSystems


"""
Whip up a DynamicalSystem from a Process
"""
function process2ds(P::Process)
    prob = process2problem(P)
    d = ContinuousDynamicalSystem(prob)
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
export DynamicalSystems.lyapunov