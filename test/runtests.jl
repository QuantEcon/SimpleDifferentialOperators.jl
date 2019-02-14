using SimpleDifferentialOperators
using Test, LinearAlgebra

@testset "Dispatch Tests" begin
    methodList = collect(methods(SimpleDifferentialOperators._diffusionoperators))

    # Test for defined cases
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractRange, Reflecting, Reflecting)) == methodList[1]
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractRange, Mixed, Mixed)) == methodList[2]
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractArray, Reflecting, Reflecting)) == methodList[3]
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractArray, Mixed, Mixed)) == methodList[4]

    # Test for error handling
    grids = [range(1.0, 10.0, length = 100), collect(range(1.0, 10.0, length = 100))]
    BCs = [Reflecting(), Mixed(2.0), Absorbing(1.0, 2.0)] # list of all possible BCs
    @test_throws MethodError diffusionoperators(grids[1], BCs[1], BCs[3])
    @test_throws MethodError diffusionoperators(grids[2], BCs[1], BCs[2])
end

@testset "Reflecting Barrier Tests" begin
    σ = 1; μ = -1;
    # Uniform grid
    L_1_minus, L_1_plus, L_2 = diffusionoperators(1:5, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(1:5, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
    @test -μ * L_1_plus + σ^2/2 * L_2 == [-1.5 1.5 0.0 0.0 0.0; 0.5 -2.0 1.5 0.0 0.0; 0.0 0.50 -2.0 1.50 0.0; 0.0 0.0 0.50 -2.0 1.50; 0.0 0.0 0.0 0.50 -0.50]
    # Irregular grid 
    L_1_minus, L_1_plus, L_2 = diffusionoperators(collect(1:5), Reflecting(), Reflecting()) # irregular grid
    @test @inferred(diffusionoperators(collect(1:5), Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
    @test -μ * L_1_plus + σ^2/2 * L_2 == [-1.5 1.5 0.0 0.0 0.0; 0.5 -2.0 1.5 0.0 0.0; 0.0 0.50 -2.0 1.50 0.0; 0.0 0.0 0.50 -2.0 1.50; 0.0 0.0 0.0 0.50 -0.50]
end
