using DynamicalSystems
using Plots
plotly()

# Simulate a sine wave using DynamicalSystems, compare accuracy and timing to builtin sine and use this to construct Process type containing (only) relevant metadata

X0 = [0.0]
T = 2000.0
μ = 1.0 # Hz

function dynamicalSine(X, μ, t)
    dX = μ.*cos(μ.*t) # dx/dt
    return SVector{1}(dX)
end
function dynamicalSineJ(X, μ, t)
    J = @SMatrix [-μ^2*sin(μ.*t)]
end

#(x, μ, t) -> sin(μ.*t)
ds = ContinuousDynamicalSystem(dynamicalSine, X0, μ, dynamicalSineJ) # Easy system defined anonymously
data = trajectory(ds, T; dt=0.01)

plot(data[:, 1], seriestype=:line)
plot!(title="Sine")