using DynamicalSystems
using NonstationaryProcessesBase
using Plots

# Integrate an FM wave with DynamicalSystems

X0 = [0.0]
T = 2000.0
μ = unitStep(1000).*1 # Hz

function fmWave(X, μ, t)
    dX = μ.*cos(μ.*t) # dx/dt
    return SVector{1}(dX)
end
# function fmWaveJ(X, μ, t)
#     J = @SMatrix [-μ^2*sin(μ.*t)]
# end

#(x, μ, t) -> sin(μ.*t)
ds = ContinuousDynamicalSystem(dynamicalSine, X0, μ, dynamicalSineJ) # Easy system defined anonymously
data = trajectory(ds, T; dt=0.01)

plot(data[:, 1], seriestype=:path)
plot!(title="Sine")
