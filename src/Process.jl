# How best to structure the I/O format and system definitions...

Base.@kwdef struct Process # Not ensemble
    process = []
    parameter_profile::Union{Function, Tuple, Array} = constantParameter # Can be a tuple of symbols, if the system has more than one parameter
    parameter_profile_parameters::Union{Tuple, Array} = [0] # Can be a tuple of tuples
    X0::Vector = [0.0, 0.0]
    t0::Union{Float64, Int64} = -10.0
    savet0::Union{Float64, Int64} = 0.0
    dt::Union{Float64, Int64} = 0.001
    savedt::Union{Float64, Int64} = 0.01
    tmax::Union{Float64, Int64} = 100.0
    alg::SciMLBase.SciMLAlgorithm = RK4()
    solver_opts::Dict = Dict(:adaptive=>false)
    parameter_rng::UInt64 = seed()
    solver_rng::UInt64 = seed()
end

# ------------------------------------------------------------------------------------------------ #
#               A function to handle simulations that are specified with a Process type            #
# ------------------------------------------------------------------------------------------------ #
simulate(P::Process) = P.process(P)
export simulate

timeseries(s::SciMLBase.AbstractTimeseriesSolution, dim::Real) = s[dim, :]
timeseries(s::SciMLBase.AbstractTimeseriesSolution, dim::Union{Vector, UnitRange}=1:size(s.u[1], 1)) = s[dim, :]'
# function timeseries(s::Tuple, dim::Union{Real, Vector, Tuple}=1)
#     timeseries(s[1], dim) # You gave the metadata as well
# end
export timeseries
