using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Operators under reflecting barrier conditions" begin
    @testset "Accuracy & regression test" begin
        # Uniform grid
        x = 1:3
        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(x, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
        @test @inferred(diffusionoperators(x, (Mixed(0.), Mixed(0.)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L₂) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -1.]
        @test x̄ == [0; x; 4]

        # Irregular grid
        x = collect(x)
        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(x, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x,(Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
        @test @inferred(diffusionoperators(x, (Mixed(0.), Mixed(0.)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L₂) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -1.]
        @test x̄ == [0; x; 4]
    end

    @testset "Consistency" begin
        # Test for consistency
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid,(Reflecting(), Reflecting()))
        L₁₋_ir, L₁₊_ir, L₂_ir, x̄_ir = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
        @test L₁₋ ≈ L₁₋_ir
        @test L₁₊ ≈ L₁₊_ir
        @test L₂ ≈ L₂_ir
        @test x̄ ≈ x̄_ir
    end
end

@testset "Operators under mixed barrier conditions" begin
    @testset "Accuracy & regression test" begin
        # Uniform grid
        x = 1:3
        ξ_lb, ξ_ub = (1., 2.)

        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(x, (Mixed(ξ_lb), Mixed(ξ_ub)))
        @test @inferred(diffusionoperators(x, (Mixed(ξ_lb), Mixed(ξ_ub)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
        @test Array(L₁₋) == [-ξ_lb 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -ξ_ub]
        @test Array(L₂) == [(-1. + ξ_lb) 1. 0.; 1. -2. 1.; 0. 1. (-1. - ξ_ub)]
        @test x̄ == [0; x; 4]

        # Irregular grid
        x = collect(x)
        ξ_lb, ξ_ub = (1., 2.)

        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(x, (Mixed(ξ_lb), Mixed(ξ_ub)))
        @test @inferred(diffusionoperators(x, (Mixed(ξ_lb), Mixed(ξ_ub)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
        @test Array(L₁₋) == [-ξ_lb 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -ξ_ub]
        @test Array(L₂) == [(-1. + ξ_lb) 1. 0.; 1. -2. 1.; 0. 1. (-1. - ξ_ub)]
        @test x̄ == [0; x; 4]
    end

    @testset "Consistency" begin
        ξ_lb, ξ_ub = (1., 2.)
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Mixed(ξ_lb), Mixed(ξ_ub)))
        L₁₋_ir, L₁₊_ir, L₂_ir, x̄_ir = diffusionoperators(irregularGrid, (Mixed(ξ_lb), Mixed(ξ_ub)))
        @test L₁₋ ≈ L₁₋_ir
        @test L₁₊ ≈ L₁₊_ir
        @test L₂ ≈ L₂_ir
        @test x̄ ≈ x̄_ir
    end
end

@testset "Operators under upwind schemes" begin
    f(x) = x^2
    μ(x) = -x # drift depends on state
    σ = 0.1
    ρ = 0.05
    M = 100 # size of grid
    x = range(-1.0, 1.0, length = M) # grid
    ## M-vector of drifts stacked according to the states
    μs = μ.(x)

    # operators with reflecting boundary conditions
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(x, (Reflecting(), Reflecting()))

    # Define first order differential operator using upwind scheme
    L_1_upwind = (μs .<= 0) .* L₁₋ + (μs .> 0) .* L₁₊

    indices_m = findall(μs .<= 0)
    indices_p = findall(μs .> 0)

    @test L_1_upwind[indices_m, indices_m] == Array(L₁₋[indices_m, indices_m])
    @test L_1_upwind[indices_p, indices_p] == Array(L₁₊[indices_p, indices_p])
end

@testset "Input Type Variance" begin
    # BigFloat
    uniformGrid = range(BigFloat(0.0), BigFloat(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(irregularGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)

    # Float32
    uniformGrid = range(Float32(0.0), Float32(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Mixed(one(Float32)), Mixed(one(Float32))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(irregularGrid, (Mixed(one(Float32)), Mixed(one(Float32))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)

    # Duals
    uniformGrid = range(Dual(0.0), Dual(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(uniformGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
    L₁₋, L₁₊, L₂, x̄ = diffusionoperators(irregularGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂, x̄ = x̄)
end
