import Base.:+
import Base.:-
import Base.:*
import Base.:/
import Base.:^

"""
     Discontinuous

 A type that contains a discontinuous function and an index of any of its discontinuities. This allows the construction of step parameter profiles that retain discontinuity information. You can also call the Discontinuous type like a normal function, if you want to.
"""
struct Discontinuous <: Function
    f::Function
    d::Set
end

arithmetics = (:+, :-, :*, :/, :^)
for a ∈ arithmetics
    eval(quote
        ($a)(c::Real, D1::Discontinuous) = Discontinuous(x -> ($a)(c, D1.f(x)), D1.d)
        ($a)(D1::Discontinuous, c::Real, ) = Discontinuous(x -> ($a)(D1.f(x), c), D1.d)
        ($a)(D1::Discontinuous, D2::Discontinuous) = Discontinuous((x...) -> ($a)(D1.f(x...), D2.f(x...)), union(D1.d, D2.d))
        ($a)(f::Function, D::Discontinuous) = Discontinuous((x...) -> ($a)(f(x...), D.f(x...)), D.d)
        ($a)(D::Discontinuous, f::Function) = Discontinuous((x...) -> ($a)(D.f(x...), f(x...)), D.d)
    end)
end


# Overload call so that you can use a Discontinuous like a normal (vectorised) function, if you want
(D::Discontinuous)(x::Real) = D.f(x)
(D::Discontinuous)(x::Union{Array, StepRange, StepRangeLen, UnitRange}) = D.f.(x)


# ! To use stiff solvers that error about instability with discontinuous profiles, try e.g. Rosenbrock23(autodiff = false) and turning off adaptive timestepping
