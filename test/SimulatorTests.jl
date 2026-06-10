using NonstationaryProcesses
using Test
using NonstationaryProcesses: timeseries, simulate  # disambiguate from Makie/other exports when CairoMakie is loaded (full suite)

# Assertion-based unit tests for the core (DifferentialEquations-free) simulators
# defined in src/Noise.jl and src/Signals.jl. These complement the @test_nowarn
# smoke scripts in runtests.jl by checking output shape, finiteness, and the
# reproducibility contract of the Process/simulate pipeline.

# Simulators whose process function calls seed(P.solver_rng): a fixed Process
# reproduces its series, while an independently constructed Process draws a fresh rng.
const REPRODUCIBLE_SIMS = (:arSim => arSim,
    :gaussianBimodalSim => gaussianBimodalSim,
    :bimodalSwitchingSim => bimodalSwitchingSim,
    :shiftyNoiseSim => shiftyNoiseSim,
    :noisySineSim => noisySineSim,
    :noisyShiftyScalySineSim => noisyShiftyScalySineSim,
    :noisyShiftyShcalySineSim => noisyShiftyShcalySineSim,
    :shcalySineSim => shcalySineSim)

# fmWave omits the seed(P.solver_rng) call, so its (stochastic stepNoise) parameter
# profile is not reproducible across simulate calls; only finiteness/variation hold.
const UNSEEDED_SIMS = (:fmWaveSim => fmWaveSim,)

@testset "reproducible simulators" begin
    for (name, sim) in REPRODUCIBLE_SIMS
        @testset "$name" begin
            P = sim()
            x = timeseries(simulate(P))
            @test !isempty(x)
            @test eltype(x) <: Real
            @test all(isfinite, x)
            @test timeseries(simulate(P)) == x          # a fixed Process reproduces its series
            @test timeseries(simulate(sim())) != x      # a fresh Process draws a new rng
        end
    end
end

@testset "unseeded simulators" begin
    for (name, sim) in UNSEEDED_SIMS
        @testset "$name" begin
            P = sim()
            x = timeseries(simulate(P))
            @test !isempty(x)
            @test eltype(x) <: Real
            @test all(isfinite, x)
            @test timeseries(simulate(P)) != x          # no seed() call -> varies between runs
        end
    end
end

@testset "tmax controls series length" begin
    for n in (100, 250, 1000)
        @test length(timeseries(simulate(arSim(tmax=n)))) == n + 1   # t0=0, savedt=1 -> tmax+1 samples
    end
end

@testset "gaussianBimodal scalar draw" begin
    # Bare (unexported) process function: a single mixture draw, always a finite real.
    for (μ, σ, α) in ((0.0, 1.0, 0.5), (5.0, 2.0, 0.2), (-3.0, 0.5, 0.8))
        x = NonstationaryProcesses.gaussianBimodal(μ, σ, α)
        @test x isa Real
        @test isfinite(x)
    end
end

# --- DifferentialEquations-based simulators (now part of the package) ---

# Continuous ODE flows that integrate deterministically: each produces a finite,
# multivariate (time x state) series, reproducible across simulate calls. Only the
# flows that simulate cleanly with a short tmax are covered here; several other
# flow/map simulators have pre-existing issues (out-of-place process functions,
# missing helpers) that are out of scope for these tests.
const FLOW_SIMS = (:lorenzSim => lorenzSim,
    :diffusionlessLorenzSim => diffusionlessLorenzSim,
    :doubleScrollSim => doubleScrollSim,
    :thomasCyclicallySymmetricSim => thomasCyclicallySymmetricSim,
    :piecewiseLinearHyperchaosSim => piecewiseLinearHyperchaosSim,
    :waveDrivenHarmonicSim => waveDrivenHarmonicSim,
    :doublePendulumSim => doublePendulumSim)

@testset "continuous flows" begin
    for (name, sim) in FLOW_SIMS
        @testset "$name" begin
            X = timeseries(simulate(sim(tmax=5.0)))
            @test !isempty(X)
            @test ndims(X) == 2                                # time x state
            @test eltype(X) <: Real
            @test all(isfinite, X)
            @test timeseries(simulate(sim(tmax=5.0))) == X     # deterministic ODE -> reproducible
        end
    end
end

@testset "Process solve infrastructure" begin
    P = lorenzSim(tmax=5.0)
    @test process2problem(P) !== nothing       # builds an ODEProblem from a Process
    sol = process2solution(P)                  # builds + solves
    @test length(sol) > 0
    @test all(isfinite, Array(sol))
end
