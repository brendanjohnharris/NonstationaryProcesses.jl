using DynamicalSystems
using Plots
using Random

# Simulate an AR process using DynamicalSystems.jl

@inline @inbounds function AR(X, p, t)
    r = randn()
    println(r)
    X = r + (p).*X
    println(X)
    return X
end
Random.seed!(32)
ds = DiscreteDynamicalSystem(AR, 0.0, 0.1)

# ----------- Look more closely at the integrator-----------
# function cut_step!(integ, N)
#     #println(integ.u[1])
#     for i in 1:N
#         integ.u = integ.f(integ.u, integ.p, integ.t)
#         integ.t += 1
#     end
#     return
# end

# function cut_trajectory(ds, t, u = ds.u0; dt = 1, Ttr = 0) where {S}
#     dt = round(Int, dt)
#     ti = ds.t0
#     tvec = ti:dt:t+ti
#     L = length(tvec)
#     integ = integrator(ds, u)
#     Ttr != 0 && step!(integ, Ttr)
#     data = zeros(L)
#     data[1] = integ.u[1]
#     for i in 2:L
#         #print(integ)
#         cut_step!(integ, dt)
#         data[i] = integ.u[1]
#     end
#     return data
# end
#data = cut_trajectory(ds, 25; dt=1)
#----------------------------------------
data = trajectory(ds, 25; dt=1)

plot(data[:, 1], seriestype=:line)

# Check that this EOM is in the correct form by reproducing:
f(X, p) = randn() + p.*X
function integrateAR(X0, p, T)
    X = zeros(T+1)
    X[1] = X0
    for t = 1:T
        X[t+1] = f(X[t], p)
    end
    return X
end
Random.seed!(32)
t = integrateAR(0.0, 0.1, 25)
plot!(t, seriestype=:line)


