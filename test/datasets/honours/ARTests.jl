using DynamicalSystems
using CairoMakie
using Random
import DynamicalSystems.trajectory
# Simulate an AR process using DynamicalSystems.jl
X0 = [0.0, 0.1]
p = [0.1, 0.1]
T = 25
@inline @inbounds function AR(X, p, t)
    r = randn()
    dX1 = r + p[1] .* X[1] + X[2]
    dX2 = p[2] .* X[1]
    return SVector{2}(dX1, dX2)
end
Random.seed!(32)
ds = DiscreteDynamicalSystem(AR, X0, p)

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
data = trajectory(ds, T)

plot(data[1][:, 1])

# Check that this EOM is in the correct form by reproducing:
function f(X, p)
    r = randn()
    dX1 = r + p[1] .* X[1] + X[2]
    dX2 = p[2] .* X[1]
    return [dX1, dX2]
end

function integrateAR(X0, p, T)
    X = zeros(T + length(p), 2)

    X[1, :] = X0
    for t = 1:size(X, 1)-1
        X[t+1, :] = f(X[t, :], p)
    end
    return X[length(p):end, :]
end
Random.seed!(32)
t = integrateAR(X0, p, T)
plot!(t[:, 1])

current_figure()
