using Distributions
# ------------------------------------------------------------------------------------------------ #
#                                         Gaussian Bimodal                                         #
# ------------------------------------------------------------------------------------------------ #

function gaussianBimodal(μ=0.0, σ=1.0, α=0.5)
    # μ is the mean of the satellite gaussian, σ is the width (SD) and α is the proportional probability (i.e. the mass of the satellite compared to the central gaussian)
    if rand() > α # Draw from the first distribution
        x = randn()
    else # Draw from the satellite
        x = σ*randn() + μ
    end
end

function gaussianBimodal(P::Process)
    # Sample from a distribution formed from a single unit gaussian at the origin, and a satellite gaussian with parameters (height, mean, width) you choose.
    seed(P.solver_rng)
    sol = [gaussianBimodal(parameter_function(P)(t)...) for t in P.transient_t0:P.savedt:P.tmax]
end


# ------------------------------------------------------------------------------------------------ #
#                                         Bimodal Switching                                        #
# ------------------------------------------------------------------------------------------------ #
function bimodalSwitching(t⃗, α::Function, δ::Function)
    v = zeros(length(t⃗))
    for t ∈ 2:length(t⃗)
        Y = rand(Binomial(1, α(t⃗[t])))
        v[t] = Y*v[t-1] + (1-Y)*(1-v[t-1])
    end
    x = randn(length(t⃗)) .+ v.*δ.(t⃗)
end

function bimodalSwitching(P::Process)
    # Draw from one Gaussian but at each time step have a probability 1-α of switching to another gaussian at distance δ away
    seed(P.solver_rng)
    p1, p2 = [P.parameter_profile[i](P.parameter_profile_parameters[i]...) for i ∈ 1:lastindex(P.parameter_profile)]
    sol = bimodalSwitching(P.transient_t0:P.savedt:P.tmax, p1, p2)
end


# ------------------------------------------------------------------------------------------------ #
#                                           Shifty Noise                                           #
# ------------------------------------------------------------------------------------------------ #
function shiftyNoise(P::Process)
    # Parameters (η, C)
    seed(P.solver_rng)
    (η, C) = parameter_functions(P)
    sol = [η(t)*randn() + C(t) for t in P.transient_t0:P.savedt:P.tmax]
end


# ------------------------------------------------------------------------------------------------ #
#                                                AR                                                #
# ------------------------------------------------------------------------------------------------ #
function AR(P::Process)
    # Variable number of AR parameters (ϕ₁, ϕ₂, ϕ₃, ...)
    seed(P.solver_rng)
    ϕ⃗ = parameter_function(P)
    p = length(ϕ⃗(gett0(P)))
    ξ⃗ = randn(p+1)
    X = zeros(length(times(P)), p+1)
    X[1, 1:length(getX0(P))] = getX0(P)
    for t ∈ 2:length(times(P))
        X[t, :], ξ⃗ = AR(X[t-1, :], ξ⃗, forcevec(ϕ⃗(t)))
    end
    return X[Int.(times(P, transient=true) .- (gett0(P)-1)), 1:length(getX0(P))]
end