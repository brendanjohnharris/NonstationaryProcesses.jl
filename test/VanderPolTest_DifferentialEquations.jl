using DifferentialEquations
using Plots
plotly()
# Plot a Van der Pol oscillator: \frac{d^2 x}{dt^2} - \mu(1-x^2) \frac{dx}{dt} + x = 0

function f!(dx, x, μ, t) # Can either be in-place f, in which case the first argument is dx, or out of place, in which case the first argument is x and the output is dx
    dx[1] = x[2]
    dx[2] = μ(t)*(1-x[1]^2).*x[2] - x[1]
end

X0 = [1.0, 1.0]
tspan = (0.0, 2000.0)
μ = unitStep(1000.0).*20.0 # A step from 0 to 20 at t = 1000
prob = ODEProblem(f!, X0, tspan, μ)
sol = solve(prob; tstops=sort(collect(μ.d)), saveat=0.001, reltol=1e-6)

plot(sol, vars=(1))
plot!(title="Van der Pol, parameter steps from 0 to 20")
