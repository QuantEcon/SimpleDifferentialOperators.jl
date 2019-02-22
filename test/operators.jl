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

@testset "Reflecting Barrier Tests" begin
    #=
        Correctness tests
    =#
    σ = 1; μ = -1;
    # Uniform grid
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(1:5, Reflecting(), Reflecting())
    @test @inferred(diffusionoperators(1:5, Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test @inferred(diffusionoperators(1:5, Mixed(0.), Mixed(0.))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar) # test that the mixed case properly nests the reflecting case
    @test μ * L_1_minus + σ^2/2 * L_2 == [-0.5 0.5 0.0 0.0 0.0; 1.5 -2.0 0.5 0.0 0.0; 0.0 1.5 -2.0 0.5 0.0; 0.0 0.0 1.5 -2.0 0.5; 0.0 0.0 0.0 1.5 -1.5]
    @test -μ * L_1_plus + σ^2/2 * L_2 == [-1.5 1.5 0.0 0.0 0.0; 0.5 -2.0 1.5 0.0 0.0; 0.0 0.50 -2.0 1.50 0.0; 0.0 0.0 0.50 -2.0 1.50; 0.0 0.0 0.0 0.50 -0.50]
    # Irregular grid
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(collect(1:5), Reflecting(), Reflecting()) # irregular grid
    @test @inferred(diffusionoperators(collect(1:5), Reflecting(), Reflecting())) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test @inferred(diffusionoperators(collect(1:5), Mixed(0.), Mixed(0.))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test μ * L_1_minus + σ^2/2 * L_2 == [-0.5 0.5 0.0 0.0 0.0; 1.5 -2.0 0.5 0.0 0.0; 0.0 1.5 -2.0 0.5 0.0; 0.0 0.0 1.5 -2.0 0.5; 0.0 0.0 0.0 1.5 -1.5]
    @test -μ * L_1_plus + σ^2/2 * L_2 == [-1.5 1.5 0.0 0.0 0.0; 0.5 -2.0 1.5 0.0 0.0; 0.0 0.50 -2.0 1.50 0.0; 0.0 0.0 0.50 -2.0 1.50; 0.0 0.0 0.0 0.50 -0.50]

    #=
        Consistency tests
    =#
    uniformGrid = range(0.0, 1.0, length = 500)
    irregularGrid = collect(uniformGrid)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Reflecting(), Reflecting())
    L_1_minus_ir, L_1_plus_ir, L_2_ir, x_bar_ir = diffusionoperators(irregularGrid, Reflecting(), Reflecting())
    @test L_1_minus ≈ L_1_minus_ir
    @test L_1_plus ≈ L_1_plus_ir
    @test L_2 ≈ L_2_ir
    @test x_bar ≈ x_bar_ir
end

@testset "Mixed Boundary Tests" begin
    σ = 1; μ = -1;
    uniformGrid = 1:1:5
    irregularGrid = collect(uniformGrid)
    ξ_1, ξ_2 = (1., 2.)
    #=
        Accuracy tests
    =#
    ξ = ξ_1
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Mixed(ξ), Mixed(ξ))
    @test @inferred(diffusionoperators(uniformGrid, Mixed(ξ), Mixed(ξ))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test μ * L_1_minus + σ^2/2 * L_2 == [-1+(1+ξ)+(-2+1+ξ)/2 0.5 0.0 0.0 0.0; 1.5 -2.0 0.5 0.0 0.0; 0.0 1.5 -2.0 0.5 0.0; 0.0 0.0 1.5 -2.0 0.5; 0.0 0.0 0.0 1.5 -1+(-2+1-ξ)/2]
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(irregularGrid, Mixed(ξ), Mixed(ξ))
    @test @inferred(diffusionoperators(irregularGrid    , Mixed(ξ), Mixed(ξ))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test μ * L_1_minus + σ^2/2 * L_2 == [-1+(1+ξ)+(-2+1+ξ)/2 0.5 0.0 0.0 0.0; 1.5 -2.0 0.5 0.0 0.0; 0.0 1.5 -2.0 0.5 0.0; 0.0 0.0 1.5 -2.0 0.5; 0.0 0.0 0.0 1.5 -1+(-2+1-ξ)/2]

    ξ = ξ_2
    L_1_minus, L_1_plus, L_2 = diffusionoperators(uniformGrid, Mixed(ξ), Mixed(ξ))
    @test @inferred(diffusionoperators(uniformGrid, Mixed(ξ), Mixed(ξ))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test μ * L_1_minus + σ^2/2 * L_2 == [-1+(1+ξ)+(-2+1+ξ)/2 0.5 0.0 0.0 0.0; 1.5 -2.0 0.5 0.0 0.0; 0.0 1.5 -2.0 0.5 0.0; 0.0 0.0 1.5 -2.0 0.5; 0.0 0.0 0.0 1.5 -1+(-2+1-ξ)/2]
    L_1_minus, L_1_plus, L_2 = diffusionoperators(irregularGrid, Mixed(ξ), Mixed(ξ))
    @test @inferred(diffusionoperators(irregularGrid, Mixed(ξ), Mixed(ξ))) == (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
    @test μ * L_1_minus + σ^2/2 * L_2 == [-1+(1+ξ)+(-2+1+ξ)/2 0.5 0.0 0.0 0.0; 1.5 -2.0 0.5 0.0 0.0; 0.0 1.5 -2.0 0.5 0.0; 0.0 0.0 1.5 -2.0 0.5; 0.0 0.0 0.0 1.5 -1+(-2+1-ξ)/2]

    #=
        Consistency tests
    =#
    uniformGrid = range(0.0, 1.0, length = 500)
    irregularGrid = collect(uniformGrid)
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(uniformGrid, Mixed(ξ_1), Mixed(ξ_2))
    L_1_minus_ir, L_1_plus_ir, L_2_ir, x_bar_ir = diffusionoperators(irregularGrid, Mixed(ξ_1), Mixed(ξ_2))
    @test L_1_minus ≈ L_1_minus_ir
    @test L_1_plus ≈ L_1_plus_ir
    @test L_2 ≈ L_2_ir
    @test x_bar ≈ x_bar_ir
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
