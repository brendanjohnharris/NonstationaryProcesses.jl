# ------------------------------------------------------------------------------------------------ #
#                                             Henon map                                            #
# ------------------------------------------------------------------------------------------------ #
@inline @inbounds function henon(x, p, n) # In this case, p is a tuple of functions
    dx1 = 1.0 - p(n)[1]*x[1]^2 + x[2]
    dx2 = p(n)[2]*x[1]
    return SVector{2}(dx1, dx2)
end
# Jacobian:
@inline @inbounds function henon_J(x, p, n)
    J = @SMatrix [-2*p(n)[1]*x[1]      1.0;
                p(n)[2]            0.0]
    return J
end

function henon(P::Process)
    seed(P.solver_rng)
    prob = ODEProblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=henon_J)
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end
