using Distributions
using PyCall
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

gaussianBimodalSim = Process(
    process = gaussianBimodal,
    parameter_profile = (ramp, constantParameter, constantParameter),
    parameter_profile_parameters = ((0.001, 0.1, 0.0), (1.0,), (0.1)),
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 1,
    tmax = 10000.0)
export gaussianBimodalSim


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

bimodalSwitchingSim = Process(
    process = bimodalSwitching,
    parameter_profile = (constantParameter, constantParameter),
    parameter_profile_parameters = ((0.9,), (10.0)),
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 1,
    tmax = 1000.0)
export bimodalSwitchingSim


# ------------------------------------------------------------------------------------------------ #
#                                           Shifty Noise                                           #
# ------------------------------------------------------------------------------------------------ #
function shiftyNoise(P::Process)
    # Parameters (η, C)
    seed(P.solver_rng)
    (η, C) = parameter_functions(P)
    sol = [η(t)*randn() + C(t) for t in P.transient_t0:P.savedt:P.tmax]
end


shiftyNoiseSim = Process(
    process = shiftyNoise,
    parameter_profile = (constantParameter, stepNoise),
    parameter_profile_parameters = [(1.0,), ((0.0, 1000.0), 100.0, 2.0, 0.0)],
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 1,
    tmax = 1000.0)
export shiftyNoiseSim


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

arSim = Process(
    process = AR,
    X0 = [0.0], # Number of initial conditions should really be 1 + num. of parameters, but if you do not specify they default to 0.0
    parameter_profile = (constant, constant, constant, constant),
    parameter_profile_parameters = ((0.2,), (0.2,), (0.2,), (0.2,)),
    t0 = 0,
    dt = 1,
    savedt = 1,
    tmax = 5000)
export arSim


"""
Wrap the main function of the python package colorednoise
"""
function colorednoise(β, N)
    cn.powerlaw_psd_gaussian(β, N)
end
export colorednoise

# """
# Colored noise with a power-law spectrum described by the exponent β
# """
# function coloredNoise(P::Process)
#     # Parameter β
#     seed(P.solver_rng) # Does this carry over to python? Not a huge deal if not
#     (β,) = parameter_functions(P)
#     sol = [η(t)*randn() + C(t) for t in P.transient_t0:P.savedt:P.tmax]
# end


shiftyNoiseSim = Process(
    process = shiftyNoise,
    parameter_profile = (constantParameter, stepNoise),
    parameter_profile_parameters = [(1.0,), ((0.0, 1000.0), 100.0, 2.0, 0.0)],
    transient_t0 = 0.0,
    t0 = 0.0,
    savedt = 1,
    tmax = 1000.0)
export shiftyNoiseSim