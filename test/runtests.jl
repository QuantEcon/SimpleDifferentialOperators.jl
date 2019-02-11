using SimpleDifferentialOperators
using Test, LinearAlgebra

@testset "Grid dispatch tests" begin
    grid = range(1., 10., length = 100) # calls the regular method
    @test which(rescaled_diffusionoperators, (typeof(grid), typeof(0.))).sig == Tuple{typeof(rescaled_diffusionoperators),AbstractRange,Any}
    grid = collect(grid) # will call the irregular method
    @test which(rescaled_diffusionoperators,  (typeof(grid), typeof(0.))).sig == Tuple{typeof(rescaled_diffusionoperators),AbstractArray,Any}
end


@testset "Operator discretization" begin
# Rescaled
    ξ_1 = 1.0
    ξ_2 = 2.0
    z_uniform = 1:5
    z_irregular = [1, 2, 3, 4, 5] # This is an AbstractArray and not a AbstractRange, so it calls the right method.

    # Test for uniform grid code, ξ_1.
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_uniform,ξ_1) # Dispatches on the discrete code.
    row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now
    # Test for irregular grid code, ξ_1.
    # Test irregular grid function produce proper regular grid
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_irregular,ξ_1) # Dispatches on the discrete code.
    row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now

    # Uniform ... ξ_2.
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_uniform,ξ_2) # Dispatches on the discrete code.
    row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now
    # Irregular ... ξ_2.
    σ = 1; μ = -1;
    z,  L_1_minus, L_1_plus, L_2 = rescaled_diffusionoperators(z_irregular,ξ_2) # Dispatches on the discrete code.
    row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 == cat(row1, row2, row3, row4, row5, dims = 1) # only test for backwards now
end
