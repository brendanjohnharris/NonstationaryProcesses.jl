import Base.:+, Base.:-, Base.:*, Base.:/

"""
     Discontinuous

 A type that contains a discontinuous function and an index of any of its discontinuities. This allows the construction of step parameter profiles that retain discontinuity information. You can also call the Discontinuous type like a normal function, if you want to.
"""
struct Discontinuous <: Function
    f::Function
    d::Set
end
# These can be turned into macros later?
function +(D1::Discontinuous, D2::Discontinuous)
    f(x) = D1.f(x) + D2.f(x) # Assume f only has one parameter. Maybe can add this as a restriction later
    d = union(D1.d, D2.d)
    D = Discontinuous(f, d)
end
function -(D1::Discontinuous, D2::Discontinuous)
    f(x) = D1.f(x) - D2.f(x)
    d = union(D1.d, D2.d)
    D = Discontinuous(f, d)
end
function *(D1::Discontinuous, D2::Discontinuous)
    f(x) = D1.f(x)*D2.f(x)
    d = union(D1.d, D2.d)
    D = Discontinuous(f, d)
end
function /(D1::Discontinuous, D2::Discontinuous)
    f(x) = D1.f(x)/D2.f(x)
    d = union(D1.d, D2.d)
    D = Discontinuous(f, d)
end

# These functions all leave the discontinuity indices unchanged
function +(D1::Discontinuous, c::Real)
    g(x) = D1.f(x) + c
    D = Discontinuous(g, D1.d)
end
function -(D1::Discontinuous, c::Real)
    g(x) = D1.f(x) - c
    D = Discontinuous(g, D1.d)
end
function *(D1::Discontinuous, c::Real)
    g(x) = D1.f(x)*c
    D = Discontinuous(g, D1.d)
end
function /(D1::Discontinuous, c::Real)
    g(x) = D1.f(x)/c
    D = Discontinuous(g, D1.d)
end

# Maybe these aren't necessary
function +(c::Real, D1::Discontinuous)
    g(x) = D1.f(x)+c
    D = Discontinuous(g, D1.d)
end
function -(c::Real, D1::Discontinuous)
    g(x) = D1.f(x)-c
    D = Discontinuous(g, D1.d)
end
function *(c::Real, D1::Discontinuous)
    g(x) = D1.f(x)*c
    D = Discontinuous(g, D1.d)
end
function /(c::Real, D1::Discontinuous)
    g(x) = D1.f(x)/c
    D = Discontinuous(g, D1.d)
end


# Overload call so that you can use a Discontinuous like a normal (vectorised) function, if you want
(D::Discontinuous)(x::Real) = D.f(x)
(D::Discontinuous)(x::Union{Array, StepRange, StepRangeLen}) = D.f.(x)


function Plots.plot(xr, D::Discontinuous, args...)
    Plots.plot(xr, D(xr), args...)
end
function Plots.plot(D::Discontinuous, args...)
    xr = extrema(collect(D.d))
    scale = abs(-(xr...))
    if scale == 0.0
        xr = collect(D.d)[1]-1:0.001:collect(D.d)[1]+1
    else
        xr = (xr[1] - 0.2*scale):(scale/1000):(xr[2] + 0.2*scale)
    end
    Plots.plot(xr, D, args...)
    plot!(xlabel="t", ylabel="p(t)", legend=nothing)
end

