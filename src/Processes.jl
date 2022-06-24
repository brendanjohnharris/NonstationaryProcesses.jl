ENV["PYTHON"]=""
using PyCall
const cn = PyNULL()
const signal = PyNULL()
run(`$(PyCall.python) -m pip install colorednoise`)
copy!(cn, pyimport("colorednoise"))
copy!(signal, pyimport_conda("scipy.signal", "scipy"))

function seed(theSeed=nothing) # Seed the rng, but return the seed. If no, nothing or NaN argument, randomly seed rng
    if isnothing(theSeed)
        theSeed = abs(Random.rand(Int64))
    end
    Random.seed!(theSeed)
    return theSeed
end

function tuplef2ftuple(f, params)
    # turn a tuple of functions into a function of tuples
    if all(isempty.(params)) # The f's are just functions on their own, no need to add parameters
        # Be warned that you can't mix these; either use all parameter functions, or all standard functions. Don't be greedy.
        if f isa Tuple
            ds = [fi isa Discontinuous ? fi.d : [] for fi in f]
            ds = reduce(∪, ds) |> Set
            pp = Discontinuous(t->[x(t) for x in f], ds)
        else
            pp = f
        end
        return pp
    elseif f isa Tuple
        ps = Vector{Function}(undef, length(f))
        for i = 1:length(f)
            ps[i] = f[i](params[i]...)
        end
        ds = [fi isa Discontinuous ? fi.d : [] for fi in ps]
        ds = reduce(∪, ds) |> Set
        p = Discontinuous(t->map((x, g) -> g(x), fill(t, length(ps)), ps), ds) # Something like that
    else
        p = f(params...)
    end
    return p
end
export tuplef2ftuple

include("Processes/ARMA.jl")
include("Processes/Noise.jl")
include("Processes/Signals.jl")
