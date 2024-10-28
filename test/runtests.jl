using TimeseriesSurrogates
using DynamicalSystems
using DifferentialEquations
using NonstationaryProcesses
using NonstationaryProcesses.NonstationaryProcessesBase
using TimeseriesTools
using Test

@testset "Simulations" begin
    @test_nowarn include("./datasets/honours/ARTests.jl")
    @test_nowarn include("./datasets/honours/HenonTest.jl")
    @test_nowarn include("./datasets/honours/lyapunovTest.jl")
    @test_nowarn include("./datasets/honours/testGaussianBimodal.jl")
end
