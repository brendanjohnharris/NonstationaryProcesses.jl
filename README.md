# NonstationaryProcesses.jl

A Julia package for simulating and analyzing non-stationary dynamical systems, providing tools to model processes with time-varying parameters.

See also [`NonstationaryProcessesBase.jl`](https://github.com/brendanjohnharris/NonstationaryProcessesBase.jl)

## Basic usage
Let's simulate a chaotic system, the Langford/Aizawa attractor, with a time-varying parameter.
```julia
using CairoMakie
using DifferentialEquations
using NonstationaryProcesses
import NonstationaryProcesses.DifferentialEquationsExt.aizawaSim
using NonstationaryProcessesBase.TimeseriesTools

tmax = 15000.0
𝛿() = t -> 3.5 .+ 0.75 * sin.(2π * t ./ tmax)

S = aizawaSim(;
    parameter_profile=([constantParameter for _ in 1:3]..., 𝛿),
    parameter_profile_parameters=(0.95, 0.7, 0.6, ()), # (𝛼, 𝛽, 𝛾, 𝛿)
    tmax)

begin
    X = Timeseries(S) |> eachcol .|> collect
    f = Figure(size=(1080, 1080), backgroundcolor=:transparent)
    ax = Axis3(f[1, 1]; aspect=(1, 1, 1), azimuth=-π / 10, elevation=π / 5)
    hidedecorations!(ax)
    hidespines!(ax)
    trajectory!(ax, X...;
        colormode=:velocity,
        linewidth=0.1, colormap=:inferno, alpha=0.25)
    f
end
```

![An example of a nonstationary time series generated by the `Langford` system](Langford.png)

## Parameter profiles
They key functionality beyond `DifferentialEquations` is a library of parameter profiles, functions that describe the time course of a parameter across time.
These include discontinuous profiles represented by a custom `Discontinuous` type, which track the time points of discontinuities and force solvers to evaluate either side of these points for accuracy.
Parameter profiles can also be provided as arbitrary functions of time.
Below are some of the pre-defined parameter profiles:

|Profile|Description|
|---|---|
|   `constantParameter`|	Constant parameter profile|
|   `heaviside`|	Heaviside step function|
|   `sigmoid`|	Sigmoid function|
|   `unitStep`|	A constant function that undergoes a step change at a specified threshold|
|   `unitBump`| A step-like "bump" (increase and then decrease) between two thresholds|
|   `compositeBump`|	A composite function combining multiple bumps|
|   `compositeStep`|	A composition of multiple step changes at specified thresholds|
|   `stepNoise`|	A composition of random bumps around a baseline|
|   `stepRandomWalk`|	A random walk defined by successive random steps|
|   `ramp`| A linear function between two points values|
|   `rampInterval`|	A ramp with saturation before and after the given interval|
|   `rampOn`|	A ramp that saturates before the interval and extrapolates after|
|   `rampOff`|	A ramp that extrapolates prior to an interval and saturates after|
|   `sineWave`|	A sine wave with a specified period, amplitude, and baseline|
|   `triangleWave`|	A triangular wave|
|   `lorentzian`|	A Lorentzian function with specified amplitude, width, center, and baseline|

## Predefined systems
Alone, this package defines some simple signals (e.g. frequency-modulated waves) and `ARMA` processes.
Loading the [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl) package will give you access to a library of deterministic flows, chaotic maps, chaotic flows, and stochastic differential equations.

## Defining new systems
New systems can be defined by specifying the equations of motion, solver glue code, and simulation parameters.
For example, a `Process` for the Lorenz attractor can be defined as follows.
First, we define the equations of motion and their gradients:
```julia
function lorenz(dX, X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    (𝜎, 𝑟, 𝑏) = p(𝑡)
    dX[1] = 𝜎 * (𝑦 - 𝑥)
    dX[2] = -𝑥 * 𝑧 + 𝑟 * 𝑥 - 𝑦
    dX[3] = 𝑥 * 𝑦 - 𝑏 * 𝑧
end
function lorenz_J(J, X::AbstractArray, p::Function, 𝑡::Real)
    (𝑥, 𝑦, 𝑧) = X
    (𝜎, 𝑟, 𝑏) = p(𝑡)
    J .= [-𝜎 𝜎 0.0;
        𝑟-𝑧 -1.0 -𝑥;
        𝑦 𝑥 -𝑏]
end
```

Next, we define some glue code that runs the `DifferentialEquations` solvers:
```julia
function lorenz(P::Process)
    seed(P.solver_rng)
    prob = odeproblem(P.process, P.X0, (P.transient_t0, P.tmax), tuplef2ftuple(P.parameter_profile, P.parameter_profile_parameters), jac=lorenz_J)
    sol = dsolve(prob, P.alg; dt=P.dt, saveat=P.savedt, P.solver_opts...)
end
```

Finally, we create a `Process` that collects the parameters of the simulation:
```julia
lorenzSim = Process(
    process=lorenz,
    X0=[0.0, -0.01, 9.0],
    parameter_profile=(constantParameter, constantParameter, constantParameter),
    parameter_profile_parameters=(10.0, 28.0, 8 / 3), # Sprott's recomendation
    transient_t0=-100.0,
    t0=0.0,
    dt=0.001,
    savedt=0.05,
    tmax=500.0,
    alg=AutoVern9(Rodas5()),
    solver_opts=Dict(:adaptive => true, :reltol => 1e-10, :abstol => 1e-10, :maxiters => 1e7))
```

To collect the time series of this simulation, we use:
```julia
NonstationaryProcesses.timeseries(lorenzSim) # A `TimeseriesTools` `AbstractTimeSeries`
```

## Simulation parameters
The `Process` syntax has some shorthand for specifying the various options for simulating a nonstationary system:

| option | description|
|---|---|
|`process`| A `Function` with a method for Process types that performs a particular simulation|
|`parameter_profile`| A tuple of `Function`s that describe how each parameter of a system evolves over time. These can be functions, or curried functions of parameters given in `:parameter_profile_parameters`|
|`parameter_profile_parameters`| A tuple of parameters (which may also be tuples) for each `:parameter_profile`|
|`X0`| The initial conditions, given as a vector equal in length to number of variables in the system|
|`t0`| The initial time of the simulation, at which the system has the state `:X0`|
|`transient_t0`| The length of time to simulate the system, from the initial conditions `:X0` at `:t0`, that is discarded when retrieving the timeseries|
|`dt`| The time step of the simulation and solver|
|`savedt`| The sampling period of the solution's timeseries, which must be a multiple of `:dt`|
|`tmax`| The final time point of the simulation. The duration of the simulation is then `:tmax` - `:transient_t0`, and the duration of the returned time series is `:tmax` - `:t0`. See times|
|`alg`| The algorithm used to solve `DifferentialEquations.jl` processes. See e.g. the list of [ODE solvers](https://diffeq.sciml.ai/stable/solvers/ode_solve/)|
|`solver_opts`| A dictionary of additional options passed to the `:alg`. This can include solver tolerances, adaptive timesteps, or the maximum number of iterations. See the common solver options|
|`solver_rng`| An integer seed for the random number generator, set to a random number by default|
|`id`| A unique integer ID for a given `Process`|
|`date`| The date at which the `Process` was created|
|`solution`| The solution of the simulation in its native format (e.g. an ODE solution)|
|`varnames`| Dummy names for the variables of a Process, defaulting to `[:x, :y, :z, ...]`|