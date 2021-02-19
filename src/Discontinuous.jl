import Base.:+, Base.:-, Base.:*, Base.:/

"""
     Discontinuous

 A type that contains a discontinuous function and an index of any of its discontinuities. This allows the construction of step parameter profiles that retain discontinuity information. You can also call the Discontinuous type like a normal function, if you want to.
"""
struct Discontinuous{T} <: Function
    f::T
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

# Overload call so that you can use a Discontinuous like a normal function, if you want
(D::Discontinuous)(x::Real) = D.f(x)
(D::Discontinuous)(x::Array) = map.(x -> D.f(x), x)
