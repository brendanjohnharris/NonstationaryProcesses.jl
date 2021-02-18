using DynamicalSystems
using Plots
plotly()

# Plot a Van der Pol oscillator: \frac{d^2 x}{dt^2} - \mu(1-x^2) \frac{dx}{dt} + x = 0
# BUT with time varying parameters

X0 = [1.0, 1.0]
T = 2000
μ(t) = heaviside(t-1000).*20
function f(X, μ, t)
    dX1 = X[2]
    dX2 = μ(t).*(1-X[1].^2).*X[2] - X[1]
    return SVector{2}(dX1, dX2)
end

ds = ContinuousDynamicalSystem(f, X0, μ)
data = trajectory(ds, T; dt=0.001)

plot(data[:, 1], seriestype=:line)
plot!(title="Van der Pol, parameter steps from 0 to 20")