using DynamicalSystems
using Plots

# Simulate an discrete dynamical system (Henon map) using DynamicalSystems.jl

@inline @inbounds function hoop(x, p, n)
    dx1 = 1.0 - p[1]*x[1]^2 + x[2]
    dx2 = p[2]*x[1]
    return SVector{2}(dx1, dx2)
end
# Jacobian:
@inline @inbounds function hoop_jac(x, p, n)
    J = @SMatrix [-2*p[1]*x[1]      1.0;
                p[2]            0.0]
    return J
end

ds = DiscreteDynamicalSystem(hoop, [0.0; 0.0], [1.4, 0.3], hoop_jac)

data = trajectory(ds, 10000; dt=1)

plot(data[:, 1], data[:, 2], seriestype=:scatter, markersize=0.5)
