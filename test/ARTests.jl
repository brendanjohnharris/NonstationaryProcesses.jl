using DynamicalSystems
using Plots

# Simulate an AR process using DynamicalSystems.jl

@inline @inbounds function hiip(dx, x, p, n)
    dx[1] = 1.0 - p[1]*x[1]^2 + x[2]
    dx[2] = p[2]*x[1]
    return dx
end
# Jacobian:
@inline @inbounds function hiip_jac(J, x, p, n)
    J[1,1] = -2*p[1]*x[1]
    J[1,2] = 1.0
    J[2,1] = p[2]
    J[2,2] = 0.0
    return J
end

ds = DiscreteDynamicalSystem(hiip, [0.0; 0.0], [1.4, 0.3], hiip_jac)

data = trajectory(ds, 10000; dt=1)

plot(data[:, 1], data[:, 2], seriestype=:scatter, markersize=0.5)