using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Dispatch Tests" begin
    methodList = collect(methods(SimpleDifferentialOperators._diffusionoperators))

    # Test for defined cases
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractRange, Reflecting, Reflecting)) == methodList[1]
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractRange, Mixed, Mixed)) == methodList[2]
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractArray, Reflecting, Reflecting)) == methodList[1]
    @test which(SimpleDifferentialOperators._diffusionoperators, (AbstractArray, Mixed, Mixed)) == methodList[2]

    # Test for error handling
    grids = [range(1.0, 10.0, length = 100), collect(range(1.0, 10.0, length = 100))]
    BCs = [Reflecting(), Mixed(2.0), Absorbing(1.0, 2.0)] # list of all possible BCs
    @test_throws MethodError diffusionoperators(grids[1], BCs[1], BCs[3])
    @test_throws MethodError diffusionoperators(grids[2], BCs[1], BCs[2])
end

@testset "Operators under reflecting barrier conditions" begin
    @testset "Accuracy & regression test" begin
        # Uniform grid
        x = 1:3
        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())
        @test @inferred(diffusionoperators(x, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
        @test @inferred(diffusionoperators(x, Mixed(0.), Mixed(0.))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar) # test that the mixed case properly nests the reflecting case
        @test Array(L_1_minus) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L_1_plus) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L_2) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -1.]
        @test x_bar == [0; x; 4]

        # Irregular grid
        x = collect(x)
        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())
        @test @inferred(diffusionoperators(x, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
        @test @inferred(diffusionoperators(x, Mixed(0.), Mixed(0.))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar) # test that the mixed case properly nests the reflecting case
        @test Array(L_1_minus) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L_1_plus) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L_2) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -1.]
        @test x_bar == [0; x; 4]
    end

    @testset "Consistency" begin
        # Test for consistency
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Reflecting(), Reflecting())
        L_1_minus_ir, L_1_plus_ir, L_2_ir, x_bar_ir = diffusionoperators(irregularGrid, Reflecting(), Reflecting())
        @test L_1_minus ≈ L_1_minus_ir
        @test L_1_plus ≈ L_1_plus_ir
        @test L_2 ≈ L_2_ir
        @test x_bar ≈ x_bar_ir
    end
end

@testset "Operators under mixed barrier conditions" begin
    @testset "Accuracy & regression test" begin
        # Uniform grid
        x = 1:3
        ξ_lb, ξ_ub = (1., 2.)

        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Mixed(ξ_lb), Mixed(ξ_ub))
        @test @inferred(diffusionoperators(x, Mixed(ξ_lb), Mixed(ξ_ub))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
        @test Array(L_1_minus) == [-ξ_lb 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L_1_plus) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -ξ_ub]
        @test Array(L_2) == [(-1. + ξ_lb) 1. 0.; 1. -2. 1.; 0. 1. (-1. - ξ_ub)]
        @test x_bar == [0; x; 4]

        # Irregular grid
        x = collect(x)
        ξ_lb, ξ_ub = (1., 2.)

        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Mixed(ξ_lb), Mixed(ξ_ub))
        @test @inferred(diffusionoperators(x, Mixed(ξ_lb), Mixed(ξ_ub))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
        @test Array(L_1_minus) == [-ξ_lb 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L_1_plus) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -ξ_ub]
        @test Array(L_2) == [(-1. + ξ_lb) 1. 0.; 1. -2. 1.; 0. 1. (-1. - ξ_ub)]
        @test x_bar == [0; x; 4]
    end

    @testset "Consistency" begin
        ξ_lb, ξ_ub = (1., 2.)
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Mixed(ξ_lb), Mixed(ξ_ub))
        L_1_minus_ir, L_1_plus_ir, L_2_ir, x_bar_ir = diffusionoperators(irregularGrid, Mixed(ξ_lb), Mixed(ξ_ub))
        @test L_1_minus ≈ L_1_minus_ir
        @test L_1_plus ≈ L_1_plus_ir
        @test L_2 ≈ L_2_ir
        @test x_bar ≈ x_bar_ir
    end
end

@testset "Input Type Variance" begin
    # BigFloat
    uniformGrid = range(BigFloat(0.0), BigFloat(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(uniformGrid, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(irregularGrid, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Mixed(one(BigFloat)), Mixed(one(BigFloat)))
    @test @inferred(diffusionoperators(uniformGrid, Mixed(one(BigFloat)), Mixed(one(BigFloat)))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Mixed(one(BigFloat)), Mixed(one(BigFloat)))
    @test @inferred(diffusionoperators(irregularGrid, Mixed(one(BigFloat)), Mixed(one(BigFloat)))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)

    # Float32
    uniformGrid = range(Float32(0.0), Float32(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(uniformGrid, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(irregularGrid, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Mixed(one(Float32)), Mixed(one(Float32)))
    @test @inferred(diffusionoperators(uniformGrid, Mixed(one(BigFloat)), Mixed(one(BigFloat)))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Mixed(one(Float32)), Mixed(one(Float32)))
    @test @inferred(diffusionoperators(irregularGrid, Mixed(one(BigFloat)), Mixed(one(BigFloat)))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)

    # Duals
    uniformGrid = range(Dual(0.0), Dual(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(uniformGrid, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(irregularGrid, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Mixed(Dual(1.0)), Mixed(Dual(1.0)))
    @test @inferred(diffusionoperators(uniformGrid, Mixed(Dual(1.0)), Mixed(Dual(1.0)))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Mixed(Dual(1.0)), Mixed(Dual(1.0)))
    @test @inferred(diffusionoperators(irregularGrid, Mixed(Dual(1.0)), Mixed(Dual(1.0)))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
end
