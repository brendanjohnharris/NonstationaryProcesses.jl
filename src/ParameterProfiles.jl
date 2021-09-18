"""Functions for constructing parameter profiles"""
constantParameter(offset::Real=0.0) = x -> offset
export constantParameter
constant = constantParameter
export constant


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
heaviside(x::AbstractArray, stepOpt::Real=1.0) = heaviside.(x, stepOpt)

sigmoid(x::Real, width=1) = 1/(1+exp(-x/width))

"""
    unitStep(threshold::Real=0.0, baseline::Real=0.0, stepHeight::Real=1.0, stepOpt::Real=1.0)

Construct a constant function of 'x' that undergoes a step of size 'stepHeight' from a 'baseline' at 'threshold'. 'stepOpt' sets the value at 0. A Discontinuous type is returned, which can be called like a normal function but also contains the x-value for which there is a discontinuity.

# Examples
```julia-repl
julia> D = unitStep();
D([-1, 0, 1])
```
"""
function unitStep(threshold::Real=0.0, baseline::Real=0.0, stepHeight::Real=1.0, stepOpt::Real=1.0)
    Discontinuous(x -> (heaviside(x-threshold, stepOpt).*stepHeight)+baseline, Set([threshold]))
end
unitStep(T::Tuple, args...) = unitStep(T[1], args...) # For compatibility with other parameter profile syntax. If T is tuple, use T[1] as the threshold
export unitStep

# * A unit step that is smooth as butter
function unitsigmoid(threshold::Real=0.0, baseline::Real=0.0, stepHeight::Real=1.0, width::Real=0.1)
    Discontinuous(x -> (sigmoid(x-threshold, width).*stepHeight)+baseline, Set([threshold]))
end
unitsigmoid(T::Tuple, args...) = unitsigmoid(T[1], args...) # For compatibility with other parameter profile syntax. If T is tuple, use T[1] as the threshold
export unitsigmoid


# Same as unitStep, but come back down
function unitBump(T::Tuple, baseline::Real=0.0, bumpHeight::Real=1.0, stepOpt::Real=1.0)
    D = unitStep(T[1], baseline, bumpHeight, stepOpt) - unitStep(T[2], 0.0, bumpHeight, stepOpt)
end
export unitBump



function compositeStep(stepThrs::Vector, baselines::Vector, stepHeights::Vector)
    D = map(unitStep, stepThrs, baselines, stepHeights)
    D = sum(D)
end
compositeStep(T::Tuple, args...) = compositeStep(T[1], args...)
export compositeStep



"""
    stepNoise()

Compose a discontinuous function of many square bumps and troughs, following a white noise process.
    'stepHeight':   the standard deviation of the step height
    'T':            a tuple, (t0, tmax)

# Examples
```julia-repl
julia> D = stepNoise();
D([-1, 0, 1])
```
"""
function stepNoise(T::Tuple, stepWidth::Real=100, stepHeight::Real=1, baseline::Real=0)
    stepIdxs = T[1]:stepWidth:T[2]-stepWidth
    steps = stepHeight.*randn(Float64, (1, length(stepIdxs)))
    ps = sum(map((x, y) -> unitBump((x, x+stepWidth), baseline, y), stepIdxs, steps))
end
export stepNoise


function stepRandomWalk(T::Tuple, stepWidth::Real=100, stepHeight::Real=1, baseline::Real=0)
    stepIdxs = T[1]:stepWidth:T[2]
    steps = stepHeight.*randn(Float64, (1, length(stepIdxs)))
    ps = sum(map((x, y) -> unitStep((x, x+stepWidth), baseline, y), stepIdxs, steps))
end
export stepRandomWalk




function ramp(gradient::Real=1, p0::Real=0, t0::Real=0)
    x -> gradient.*(x.-t0) .+ p0
end
function ramp(p1::Real, p2::Real, t1::Real, t2::Real)
    gradient = (p2 - p1)/(t2 - t1)
    ramp(gradient, p1, t1)
end
export ramp

"""
Define a line over an interval. Saturate before and after this interval
"""
function rampInterval(p1::Real, p2::Real, t1::Real, t2::Real)
    r = ramp(p1, p2, t1, t2)
    b = unitBump((t1, t2), 0, 1)
    s1 = unitStep(t1, p1, -p1)
    s2 = unitStep(t2, 0, p2)
    return r*b + s1 + s2
end
export rampInterval

"""
Define a line over an interval. Before this interval, saturate, but after, extrapolate.
"""
function rampOn(p1::Real, p2::Real, t1::Real, t2::Real)
    r = ramp(p1, p2, t1, t2)
    b = unitStep(t1, 0, 1)
    s1 = unitStep(t1, p1, -p1)
    return r*b + s1
end
export rampOn

"""
Define a line over an interval. After this interval, saturate, but before, extrapolate.
"""
function rampOff(p1::Real, p2::Real, t1::Real, t2::Real)
    r = ramp(p1, p2, t1, t2)
    b = unitStep(t2, 1, -1)
    s2 = unitStep(t2, 0, p2)
    return r*b + s2
end
export rampOff


function sineWave(period::Real=1, amplitude::Real=1, t0::Real=0, baseline::Real=5*amplitude)
    x -> amplitude.*sin.((2π/period).*(x.-t0)) + baseline
end
export sineWave

function triangleWave(period::Real=1, amplitude::Real=1, t0::Real=0, tmax=100, baseline::Real=5*amplitude)
    # This is a discontinuous wave, so we can't (yet) define it over all x
    # Be careful to set the t0 and tmax to match the simulation time
    f = x -> (2amplitude/π).*asin.(sin.((2π/period).*(x.-t0))) + baseline
    d = Set((t0 + π/2):π:tmax) # The set of maxima and minima of the triangle wave
    Discontinuous(f, d)
end
export triangleWave



function lorentzian(A=1.0, Γ=1.0, x₀=0.0, y₀=0.0)
    # Γ is F.W.H.M, x₀ is centre, A is height and y₀ is the baseline
    # (function is not mass normalised, as would be usual)
    γ = Γ/2
    x -> (A.*γ.^2)./((x.-x₀).^2 .+ γ.^2) + y₀
end
export lorentzian
