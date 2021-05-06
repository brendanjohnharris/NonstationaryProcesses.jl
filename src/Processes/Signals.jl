# ------------------------------------------------------------------------------------------------ #
#                                           FM Signal                                              #
# ------------------------------------------------------------------------------------------------ #
function fmWave(P::Process)
    # Only parameter is the FM signal
    T = P.transient_t0:P.savedt:P.tmax
    p = parameter_functions(P)
    Δ𝑓 = 1.0
    sol = zeros(size(T))
    pint = 0.0
    for i ∈ 2:lastindex(T)
        t = T[i]
        pint += (p(t-P.savedt) + p(t))*P.savedt/2 # Crude integration, should be fine
        sol[i] = cos(2π*(t + Δ𝑓*pint))
    end
    return sol
end



# ------------------------------------------------------------------------------------------------ #
#                                            Noisy Sine                                            #
# ------------------------------------------------------------------------------------------------ #
function noisySine(P::Process)
    seed(P.solver_rng)
    sol = [sin(t) + parameter_function(P)(t)*randn() for t in P.transient_t0:P.savedt:P.tmax]
end




# ------------------------------------------------------------------------------------------------ #
#                                      Noisy Trendy Scaly Sine                                     #
# ------------------------------------------------------------------------------------------------ #
function noisyShiftyScalySine(P::Process)
    # Parameters (η, C, A)
    seed(P.solver_rng)
    (η, C, A) = parameter_functions(P)
    sol = [A(t)*sin(t) + η(t)*randn() + C(t) for t in P.transient_t0:P.savedt:P.tmax]
end

function shcalySine(P::Process)
    seed(P.solver_rng)
    A = parameter_function(P)
    sol = [A(t)*(sin(t) + 1.0*randn() + 1.0) for t in P.transient_t0:P.savedt:P.tmax]
end

