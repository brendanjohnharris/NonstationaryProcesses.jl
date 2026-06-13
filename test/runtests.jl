using TimeseriesSurrogates
using NonstationaryProcesses
using NonstationaryProcesses.NonstationaryProcessesBase
using TimeseriesBase
using Test

@testset "Simulations" begin
    @test_nowarn include("./datasets/honours/testGaussianBimodal.jl")
end

# Reference tests that cross-check against DynamicalSystems.jl. The modern
# DynamicalSystems stack only resolves on Julia >= 1.11; on 1.10 an older
# DynamicalSystems is selected, which pins StateSpaceSets to a range that
# conflicts with the SciML versions this package now depends on. To keep the
# base test environment resolvable on 1.10, DynamicalSystems and
# DifferentialEquations are omitted from test/Project.toml and added on demand
# only where the dependency graph is satisfiable.
if VERSION >= v"1.11"
    using Pkg
    Pkg.add(["DynamicalSystems", "DifferentialEquations"])
    @testset "DynamicalSystems references" begin
        @test_nowarn include("./datasets/honours/ARTests.jl")
        @test_nowarn include("./datasets/honours/HenonTest.jl")
        @test_nowarn include("./datasets/honours/lyapunovTest.jl")
    end
else
    @info "Skipping DynamicalSystems-based reference tests on Julia $VERSION (<1.11); the DynamicalSystems dependency graph is unsatisfiable there."
end

include("./SimulatorTests.jl")
