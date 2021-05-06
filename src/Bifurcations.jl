using .BifurcationKit
using .Setfield


# ------------------------------------------------------------------------------------------------ #
#                           Generate a bifurcation diagram from a Process                          #
# ------------------------------------------------------------------------------------------------ #
"""Ok, so this one plots a bifurcation diagram using the equations defining a Process.
You can give the single parameter you want to vary, param, and the range over which you want to vary it, prange. If you don't supply prange, I'll take the parameter range as the extrema of the parameter function in the Process. Note that this will only work for Processes defined using ODEproblems from DifferentialEquations. We'll discard the jacobian and compute it automatically.
"""
function bifurcation_diagram(P::Process; param=1, prange=extrema(parameterseries(P, p=param)), kwargs...)
    f = constantine(P; param)
    ode = ODEProblem(f, P.X0, (P.transient_t0, P.tmax); atol = 1e-10, rtol = 1e-9)
    opts = ContinuationPar(pMin = min(prange...),
                            pMax = max(prange...),
                            ds = 0.01, dsmax = 0.05,
	                        nInversion = 16,
                            detectBifurcation = 3,
                            maxBisectionSteps = 50,
                            nev = 3,
                            kwargs...)
    br, = continuation((x, p) -> Array(f(x, p, 0.0)), P.X0, [parameterseries(P, p=param)[1]],
                            (@lens _[1]), opts;
                            plot = true,
                            verbosity = 1)
end
export bifurcation_diagram


function constantine(P::Process; param)
    # Want to turn our special non-stationary EOM to a normal, constant parameter one
    eom = P.process
    p_profiles = forcevec(deepcopy(P.parameter_profile))
    p_parameters= forcevec(deepcopy(P.parameter_profile_parameters))
    p_profiles′ = (p_profiles[1:param-1]..., constantParameter, p_profiles[param+1:end]...)
    function f(X::AbstractArray, 𝑝, 𝑡)
        p_parameters′ = (p_parameters[1:param-1]..., 𝑝[1], p_parameters[param+1:end]...)
        if typeof(p_profiles′) <: Union{Tuple, Vector} && length(p_profiles′) == 1
            p_profiles′ = p_profiles′[1]
        end
        pf = tuplef2ftuple(p_profiles′, p_parameters′)
        eom(X, pf, 𝑡)
    end
    return f
end