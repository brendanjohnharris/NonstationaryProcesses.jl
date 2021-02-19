include("Discontinuous.jl")

# ------------------------------------------------------------------------------------------------ #
#                           Functions for constructing parameter profiles                          #
# ------------------------------------------------------------------------------------------------ #

function heaviside(x::Real, stepOpt::Real=1.0)
    if x < 0
        y = 0
    elseif x > 0
        y = 1
    elseif x == 0
        y = stepOpt
    else
        return nothing
    end
    y = convert(typeof(x), y)
end
heaviside(x::AbstractArray, stepOpt::Real=1.0) = map.(x -> heaviside(x, stepOpt), x)

"""
    unitStep(threshold::Real=0.0 stepOpt::Real=1.0)

Construct a constant function of 'x' that undergoes a unit step at 'threshold' (the Heaviside step function translated in x by 'threshold'). 'stepOpt' sets the value at 0. A Discontinuous type is returned, which can be called like a normal function but also contains the x-value for which there is a discontinuity.

# Examples
```julia-repl
julia> D = unitStep(0.0, 1.0);
D([-1, 0, 1])
```
"""

function unitStep(threshold::Real=0.0, stepOpt::Real=1.0)
    Discontinuous(x -> heaviside(x-threshold, stepOpt), Set([threshold]))
end


