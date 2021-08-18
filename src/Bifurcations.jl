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
    opts = ContinuationPar(pMin = min(prange...)|>Float64,
                            pMax = max(prange...)|>Float64,
                            ds = 0.01,
                            dsmax = 0.05,
	                        # nInversion = 16,
                            # detectBifurcation = 3,
                            # maxBisectionSteps = 50,
                            nev = 3,
                            kwargs...) # parameterseries(P, p=param)[1]
    br, = continuation((x, p) -> Array(f(x, p, 0.0)), P.X0, [min(prange...)],
                            (@lens _[1]), opts;
                            plot = true,
                            verbosity = 0)
end
export bifurcation_diagram

# jet = BifurcationKit.get3Jet(f)
# opts = ContinuationPar(pMin = min(prange...)|>Float64,
#                         pMax = max(prange...)|>Float64,
#                         ds = 0.01,
#                         dsmax = 0.05,
#                         # nInversion = 16,
#                         # detectBifurcation = 3,
#                         # maxBisectionSteps = 50,
#                         nev = 3)
# bifurcationdiagram(jet...,
#     P.X0, 10.0, (@lens _[1]),
# 	# important argument: this is the maximal
# 	# recursion level
# 	5,
#     opts)


function constantine(P::Process; param)
    # Want to turn our special non-stationary EOM to a normal, constant parameter one
    eom = P.process
    p_profiles = forcevec(deepcopy(P.parameter_profile))
    p_parameters= forcevec(deepcopy(P.parameter_profile_parameters))
    p_profiles′ = (p_profiles[1:param-1]..., constantParameter, p_profiles[param+1:end]...)
    function f(X::AbstractArray, 𝑝, 𝑡)
        p_parameters′ = (p_parameters[1:param-1]..., 𝑝[1], p_parameters[param+1:end]...)
        if p_profiles′ isa Union{Tuple, Vector} && length(p_profiles′) == 1
            p_profiles′ = p_profiles′[1]
        end
        pf = tuplef2ftuple(p_profiles′, p_parameters′)
        eom(X, pf, 𝑡)
    end
    return f
end



# optnew = NewtonPar(verbose = true, tol = 1e-12)
# opts = ContinuationPar(dsmin = 0.0001, dsmax = 0.01, ds = 0.01, pMax = 1.,
# 	newtonOptions = setproperties(optnew; maxIter = 30, tol = 1e-8),
# 	maxSteps = 300, plotEveryStep = 40,
# 	detectBifurcation = 3, nInversion = 4, tolBisectionEigenvalue = 1e-17, dsminBisection = 1e-7)

#     function optrec(x, p, l; opt = opts)
#         return setproperties(opt; maxSteps = 300, detectBifurcation = 3, detectLoop = false)
#     end

# g = (x, p) -> collect(f(x, p, 0.0))
# J = (x, p) -> x0 -> BifurcationKit.ForwardDiff.jacobian(y->g(y, p), x0)
# jet = BifurcationKit.get3Jet(g, J)
# diagram = @time bifurcationdiagram(jet...,
# 	[0.1], 1.0, (@lens _[1]),
# 	# here we specify a maximum branching level of 4
# 	4, optrec)
