using DynamicalSystems
using Plots
plotly()

# Plot a Van der Pol oscillator: \frac{d^2 x}{dt^2} - \mu(1-x^2) \frac{dx}{dt} + x = 0
# BUT with time varying parameters

X0 = [1.0, 1.0]
T = 2000.0

# --------------------------------- Generate a parameter profile --------------------------------- #
μ = unitStep(1000.0).*20.0 # A step from 0 to 20 at t = 1000
# ------------------------------------------------------------------------------------------------ #

# ----------------------------- Define an EOM, and possibly a Jacobian --------------------------- #
function vanderpol(X, μ, t)
    dX1 = X[2]
    dX2 = μ(t).*(1.0-X[1].^2).*X[2] - X[1]
    return SVector{2}(dX1, dX2)
end
# ------------------------------------------------------------------------------------------------ #

# ----------------------------------- Create a dynamical system ---------------------------------- #
ds = ContinuousDynamicalSystem(vanderpol, X0, μ)
# ------------------------------------------------------------------------------------------------ #

# --------------------------------------------- Solve -------------------------------------------- #
data = trajectory(ds, T; dt=0.001, tstops=collect(μ.d))
# ------------------------------------------------------------------------------------------------ #

plot(data[:, 1], seriestype=:line)
plot!(title="Van der Pol, parameter steps from 0 to 20")


