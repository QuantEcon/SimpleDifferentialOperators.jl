using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Operators under reflecting barrier conditions" begin

    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋ = L₁₋(x̄, bc), L₁₊ = L₁₊(x̄, bc), L₂ = L₂(x̄, bc))

    @testset "Accuracy & regression test" begin
        # Uniform grid
        x̄ = 0:4
        L₁₋, L₁₊, L₂ = diffusionoperators(x̄, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x̄, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂)
        @test @inferred(diffusionoperators(x̄, (Mixed(0.), Mixed(0.)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L₂) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -1.]

        # Irregular grid
        x̄ = [0.0; 1.0; 2.0; 4.0; 6.0]
        L₁₋, L₁₊, L₂ = diffusionoperators(x̄, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x̄,(Reflecting(), Reflecting()))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂)
        @test @inferred(diffusionoperators(x̄, (Mixed(0.), Mixed(0.)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋) == [0. 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. 0.]
        @test Array(L₂) == [-1. 1. 0.; 2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 0. 1/4 -1/4]
    end

    @testset "Consistency" begin
        # Test for consistency
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋, L₁₊, L₂ = diffusionoperators(uniformGrid,(Reflecting(), Reflecting()))
        L₁₋_ir, L₁₊_ir, L₂_ir = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
        @test L₁₋ ≈ L₁₋_ir
        @test L₁₊ ≈ L₁₊_ir
        @test L₂ ≈ L₂_ir
    end
end

@testset "Operators under mixed barrier conditions" begin

    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋ = L₁₋(x̄, bc), L₁₊ = L₁₊(x̄, bc), L₂ = L₂(x̄, bc))

    @testset "Accuracy & regression test" begin
        # Uniform grid
        x̄ = 0:4
        ξ_lb, ξ_ub = (1., 2.)

        L₁₋, L₁₊, L₂ = diffusionoperators(x̄, (Mixed(ξ_lb), Mixed(ξ_ub)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ_lb), Mixed(ξ_ub)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂)
        @test Array(L₁₋) == [-ξ_lb 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -ξ_ub]
        @test Array(L₂) == [(-1. + ξ_lb) 1. 0.; 1. -2. 1.; 0. 1. (-1. - ξ_ub)]

        # Irregular grid
        x̄ = [0.0; 1.0; 2.0; 4.0; 6.0]
        ξ_lb, ξ_ub = (1., 2.)

        L₁₋, L₁₊, L₂ = diffusionoperators(x̄, (Mixed(ξ_lb), Mixed(ξ_ub)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ_lb), Mixed(ξ_ub)))) == (L₁₋ = L₁₋, L₁₊ = L₁₊, L₂ = L₂)
        @test Array(L₁₋) == [-ξ_lb 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. -ξ_ub]
        @test Array(L₂) == [(-1. + ξ_lb) 1. 0.; 2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 0. 1/4 (-1. - ξ_ub*2.0)/4]
    end

    @testset "Consistency" begin
        ξ_lb, ξ_ub = (1., 2.)
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋, L₁₊, L₂ = diffusionoperators(uniformGrid, (Mixed(ξ_lb), Mixed(ξ_ub)))
        L₁₋_ir, L₁₊_ir, L₂_ir = diffusionoperators(irregularGrid, (Mixed(ξ_lb), Mixed(ξ_ub)))
        @test L₁₋ ≈ L₁₋_ir
        @test L₁₊ ≈ L₁₊_ir
        @test L₂ ≈ L₂_ir
    end
end

@testset "Operators under upwind schemes" begin
    f(x) = x^2
    μ(x) = -x # drift depends on state
    σ = 0.1
    ρ = 0.05
    M = 100 # size of grid
    x̄ = range(-1.0, 1.0, length = M+2) # grid
    x = x̄[2:end-1]
    bc = (Reflecting(), Reflecting()) # specify BC (reflecting barrier)
    ## M-vector of drifts stacked according to the states
    μs = μ.(x)

    # Define first order differential operator using upwind scheme
    L_1_upwind = (μs .<= 0) .* L₁₋(x̄, bc) + (μs .> 0) .* L₁₊(x̄, bc)

    indices_m = findall(μs .<= 0)
    indices_p = findall(μs .> 0)

    @test L_1_upwind[indices_m, indices_m] == Array(L₁₋(x̄, bc)[indices_m, indices_m])
    @test L_1_upwind[indices_p, indices_p] == Array(L₁₊(x̄, bc)[indices_p, indices_p])
end

@testset "Input Type Variance" begin
    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋ = L₁₋(x̄, bc), L₁₊ = L₁₊(x̄, bc), L₂ = L₂(x̄, bc))

    # BigFloat
    uniformGrid = range(BigFloat(0.0), BigFloat(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)

    # Float32
    uniformGrid = range(Float32(0.0), Float32(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Mixed(one(Float32)), Mixed(one(Float32))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋_cache, L₁₊ = L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Mixed(one(Float32)), Mixed(one(Float32))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(one(BigFloat)), Mixed(one(BigFloat))))) == (L₁₋ = L₁₋_cache, L₁₊ = L₁₊_cache, L₂ = L₂_cache)

    # Duals
    uniformGrid = range(Dual(0.0), Dual(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
    L₁₋_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(Dual(1.0)), Mixed(Dual(1.0))))) == (L₁₋ = L₁₋_cache, L₁₊ =  L₁₊_cache, L₂ = L₂_cache)
end
