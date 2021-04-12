# ------------------------------------------------------------------------------------------------ #
#                                   Skewed Pitchfork Bifurcation                                   #
# ------------------------------------------------------------------------------------------------ #

@inline @inbounds function skewedGaussianQuadratic(x::AbstractArray, p::Function, t::Real)
    (μ, η, κ) = p(t)
    x = μ*x[1] + μ*x[1]*exp(-κ*x[1])*(κ*x[1]-2.0)
    return SVector{1}(x)
end

@inline @inbounds function skewedGaussianQuadratic_σ(σ::AbstractArray, p::Function, t::Real)
    (μ, η, κ) = p(t)
    σ = η
end

function skewedGaussianQuadratic(P::Process)
    seed(P.solver_rng)
    prob = SDEProblem(P.process, skewedGaussianQuadratic_σ, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters))
    sol = dsolve(prob, P.alg; dt = P.dt, saveat=P.savedt, P.solver_opts...)
end

