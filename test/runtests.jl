using SimpleDifferentialOperators
using Test, LinearAlgebra

@testset "Discretization Tests: Consistency" begin
    ξ_vec = [1., 2., 3., 4., 5.]
    μ = -1.
    σ = 1.
    for ξ in ξ_vec
        grid = range(1., 10., length = 500)
        res1 = robin_diffusionoperators(grid, ξ)
        res2 = robin_diffusionoperators(collect(grid), ξ) # this is an AbstractArray and not a AbstractRange, so it calls the right method.
        @test μ*res1.L_1_minus + σ^2 / 2 * res1.L_2 ≈  μ*res2.L_1_minus + σ^2 / 2 * res2.L_2
        @test -μ*res1.L_1_plus + σ^2 / 2 * res1.L_2 ≈  -μ*res2.L_1_plus + σ^2 / 2 * res2.L_2
    end
end

@testset "Discretization Tests: Accuracy" begin
    ξ_1, ξ_2 = (1., 2.)
    z_uniform = 1:1:5
    z_irregular = collect(z_uniform) # dispatches properly for the same reason as L9

    # ξ_1, uniform grid.
    σ = 1; μ = -1;
    L_1_minus, L_1_plus, L_2 = robin_diffusionoperators(z_uniform, ξ_1)
    row1 = [-1 + (1+ξ_1) + (-2+1+ξ_1)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_1)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 ≈ cat(row1, row2, row3, row4, row5, dims = 1)

    # ξ_1, irregular grid.
    L_1_minus, L_1_plus, L_2 = robin_diffusionoperators(z_irregular, ξ_1)
    @test μ * L_1_minus + σ^2/2 * L_2 ≈ cat(row1, row2, row3, row4, row5, dims = 1)

    # ξ_2, regular grid.
    σ = 1; μ = -1;
    L_1_minus, L_1_plus, L_2 = robin_diffusionoperators(z_uniform,ξ_2)
    row1 = [-1 + (1+ξ_2) + (-2+1+ξ_2)/2 , 0.5, 0.0, 0.0, 0.0]'
    row2 = [1.5, -2.0, 0.5, 0.0, 0.0]'
    row3 = [0.0, 1.5, -2.0, 0.5, 0.0]'
    row4 = [0.0, 0.0, 1.5, -2.0, 0.5]'
    row5 = [0.0, 0.0, 0.0, 1.5, -1 + (-2 + 1-ξ_2)/2]'
    @test μ * L_1_minus + σ^2/2 * L_2 ≈ cat(row1, row2, row3, row4, row5, dims = 1)

    # ξ_1, irregular grid.
    L_1_minus, L_1_plus, L_2 = robin_diffusionoperators(z_irregular, ξ_2)
    @test μ * L_1_minus + σ^2/2 * L_2 ≈ cat(row1, row2, row3, row4, row5, dims = 1)
end

@testset "Tests with Varying Types" begin
    T = [Float64, BigFloat, Float32, Real] # can add DualNumber, but requires package
    for type in T
        ξ = one(type)
        grid = range(zero(type), one(type), length = 400)
        L_1_minus, L_1_plus, L_2 = robin_diffusionoperators(grid, ξ)
        @test 1 == 1  # returns true if the above call returns successfully
        @inferred robin_diffusionoperators(grid, ξ)
        @test 1 == 1 # returns true if the above inference passes
    end
end
