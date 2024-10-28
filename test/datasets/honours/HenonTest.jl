using DynamicalSystems
using CairoMakie

# Simulate an discrete dynamical system (Henon map) using DynamicalSystems.jl

@inline @inbounds function hoop(x, p, n)
    dx1 = 1.0 - p[1] * x[1]^2 + x[2]
    dx2 = p[2] * x[1]
    return SVector{2}(dx1, dx2)
end

ds = DiscreteDynamicalSystem(hoop, [0.0; 0.0], [1.4, 0.3])

data = trajectory(ds, 10000)

plot(data[1][:, 1], data[1][:, 2], markersize=1)
