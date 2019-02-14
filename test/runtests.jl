using SimpleDifferentialOperators
using Test, LinearAlgebra

@testset "Dispatch Tests" begin
    grids = [range(1.0, 10.0, length = 100), collect(range(1.0, 10.0, length = 100))]
    BCs = [Reflecting(), Mixed(2.0), Absorbing(1.0, 2.0)] # list of all possible BCs

    # Test for defined cases
    @test diffusionoperators(grids[1], BCs[1], BCs[1]) == (x = "Uniform", BC1 = "Reflecting", BC2 = "Reflecting")
    @test diffusionoperators(grids[1], BCs[2], BCs[2]) == (x = "Uniform", BC1 = "Mixed", BC2 = "Mixed")
    @test diffusionoperators(grids[2], BCs[1], BCs[1]) == (x = "Irregular", BC1 = "Reflecting", BC2 = "Reflecting")
    @test diffusionoperators(grids[2], BCs[2], BCs[2]) == (x = "Irregular", BC1 = "Mixed", BC2 = "Mixed")

    # Test for error handling
    @test_throws MethodError diffusionoperators(grids[1], BCs[1], BCs[3])
    @test_throws MethodError diffusionoperators(grids[2], BCs[1], BCs[2])
end
