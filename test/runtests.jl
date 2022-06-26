using NonstationaryProcesses
using Test

@testset "Simulations" begin
    @test_nowarn include("./ARTests.jl")
    @test_nowarn include("./HenonTest.jl")
    @test_nowarn include("./lyapunovTest.jl")
    @test_nowarn include("./testGaussianBimodal.jl")
end
