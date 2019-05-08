using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Operators under reflecting barrier conditions" begin

    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋bc = L₁₋bc(x̄, bc), L₁₊bc = L₁₊bc(x̄, bc), L₂bc = L₂bc(x̄, bc), x = interiornodes(x̄, bc))

    @testset "Accuracy & regression test" begin
        # Uniform grid
        x̄ = 0:4
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x̄, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = 0.), Mixed(ξ = 0.)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋bc) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L₂bc) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -1.]
        @test Array(x) == [1; 2; 3]

        # Irregular grid (with Δ_1m = Δ_1p and Δ_Mp = Δ_Mm)
        x̄ = [0.0; 1.0; 2.0; 4.0; 6.0]
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x̄,(Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = 0.), Mixed(ξ = 0.)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋bc) == [0. 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. 0.]
        @test Array(L₂bc) == [-1. 1. 0.; 2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 0. 1/4 -1/4]
        @test Array(x) == [1.0; 2.0; 4.0]

        # Irregular grid (with Δ_1m != Δ_1p and Δ_Mp != Δ_Mm)
        x̄ = [-1.0; 1.0; 2.0; 4.0; 7.0]

        Δ_1m = x̄[2] - x̄[1]
        Δ_1p = x̄[3] - x̄[2]
        Δ_Mm = x̄[end-1] - x̄[end-2]
        Δ_Mp = x̄[end] - x̄[end-1]
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Reflecting(), Reflecting()))
        @test @inferred(diffusionoperators(x̄,(Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test @inferred(diffusionoperators(x̄,(Mixed(ξ = 0.), Mixed(ξ = 0.)))) == diffusionoperators(x̄,(Reflecting(), Reflecting())) # test that the mixed case properly nests the reflecting case
        @test Array(L₁₋bc) == [0. 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. 0.]
        Ξ_1 = -2*(1/(Δ_1m*Δ_1p)-1/((Δ_1p+Δ_1m)*(Δ_1m)))
        Ξ_M = -2*(1/(Δ_Mm*Δ_Mp)-1/((Δ_Mp+Δ_Mm)*(Δ_Mp)))
       
        @test Array(L₂bc) == [Ξ_1 2/(Δ_1p*(Δ_1m+Δ_1p)) 0.; 
                            2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 
                            0. 2/(Δ_Mm*(Δ_Mm+Δ_Mp)) Ξ_M]
        @test Array(x) == [1.0; 2.0; 4.0]
    end

    @testset "Consistency" begin
        # Test for consistency
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(uniformGrid,(Reflecting(), Reflecting()))
        L₁₋bc_ir, L₁₊_ir, L₂_ir, x_ir = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
        @test L₁₋bc ≈ L₁₋bc_ir
        @test L₁₊bc ≈ L₁₊_ir
        @test L₂bc ≈ L₂_ir
    end
end

@testset "Operators under mixed barrier conditions" begin

    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋bc = L₁₋bc(x̄, bc), L₁₊bc = L₁₊bc(x̄, bc), L₂bc = L₂bc(x̄, bc), x = interiornodes(x̄, bc))

    @testset "Accuracy & regression test" begin
        # Uniform grid
        x̄ = 0:4
        ξ_lb, ξ_ub = (-1., -2.)

        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test Array(L₁₋bc) == [1+1/(-1+ξ_lb) 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -1+1/(1+ξ_ub)]
        @test Array(L₂bc) == [-2-1/(-1+ξ_lb) 1. 0.; 1. -2. 1.; 0. 1. -2+1/(1+ξ_ub)]
        @test Array(x) == [1; 2; 3]

        # Irregular grid (with Δ_1m = Δ_1p and Δ_Mp = Δ_Mm)
        x̄ = [0.0; 1.0; 2.0; 4.0; 6.0]
        ξ_lb, ξ_ub = (-1., -2.)
        Δ_Mm = x̄[end-1] - x̄[end-2]
        Δ_Mp = x̄[end] - x̄[end-1]

        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test Array(L₁₋bc) == [1+1/(-1+ξ_lb) 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. (-1+1/(1+ξ_ub*Δ_Mp))/Δ_Mp]
        @test Array(L₂bc) == [-1.5 1. 0.; 2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 0. 1/4 -2*(1/(Δ_Mm*Δ_Mp)-1/((1+ξ_ub*Δ_Mp)*(Δ_Mp+Δ_Mm)*(Δ_Mp)))]
        @test Array(x) == [1.0; 2.0; 4.0]

        # Irregular grid (with Δ_1m != Δ_1p and Δ_Mp != Δ_Mm)
        x̄ = [-1.0; 1.0; 2.0; 4.0; 7.0]
        ξ_lb, ξ_ub = (-1., -2.)
        Δ_1m = x̄[2] - x̄[1]
        Δ_1p = x̄[3] - x̄[2]
        Δ_Mm = x̄[end-1] - x̄[end-2]
        Δ_Mp = x̄[end] - x̄[end-1]

        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test Array(L₁₋bc) ≈ [1/3 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) ≈ [-1. 1. 0.; 0. -1/2 1/2; 0. 0. -0.4]
        @test Array(L₂bc) ≈ [-8/9 2/3 0.; 2/3 -1. 1/3; 0. 0.2 -0.36]
        @test Array(x) ≈ [1.0; 2.0; 4.0]
    end

    @testset "Accuracy & regression test, for different directions" begin
        # Uniform grid
        x̄ = 0:4
        ξ_lb, ξ_ub = (-1., -2.)

        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub, direction = :forward)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub, direction = :forward)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test Array(L₁₋bc) == [1.0 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 2.0]
        @test Array(L₂bc) == [-2.0 1. 0.; 1. -2. 1.; 0. 1. 1.0]
        @test Array(x) == [1; 2; 3]

        # Irregular grid (with Δ_1m = Δ_1p and Δ_Mp = Δ_Mm)
        x̄ = [0.0; 1.0; 2.0; 4.0; 6.0]
        ξ_lb, ξ_ub = (-1., -2.)
        Δ_Mm = x̄[end-1] - x̄[end-2]
        Δ_Mp = x̄[end] - x̄[end-1]

        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub, direction = :forward)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub, direction = :forward)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test Array(L₁₋bc) == [1.0 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. 4.0]
        @test Array(L₂bc) == [-2.0 1. 0.; 2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 0. 1/4 3/4]
        @test Array(x) == [1.0; 2.0; 4.0]

        # Irregular grid (with Δ_1m != Δ_1p and Δ_Mp != Δ_Mm)
        x̄ = [-1.0; 1.0; 2.0; 4.0; 7.0]
        ξ_lb, ξ_ub = (-1., -2.)
        Δ_1m = x̄[2] - x̄[1]
        Δ_1p = x̄[3] - x̄[2]
        Δ_Mm = x̄[end-1] - x̄[end-2]
        Δ_Mp = x̄[end] - x̄[end-1]

        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub, direction = :forward)))
        @test @inferred(diffusionoperators(x̄, (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub, direction = :forward)))) == (L₁₋bc = L₁₋bc, L₁₊bc = L₁₊bc, L₂bc = L₂bc, x = x)
        @test Array(L₁₋bc) ≈ [2.0 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) ≈ [-1. 1. 0.; 0. -1/2 1/2; 0. 0. 6.0]
        @test Array(L₂bc) ≈ [-4/3 2/3 0.; 2/3 -1. 1/3; 0. 0.2 0.6]
        @test Array(x) ≈ [1.0; 2.0; 4.0]
    end

    @testset "Consistency" begin
        ξ_lb, ξ_ub = (1., 2.)
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(uniformGrid, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))
        L₁₋bc_ir, L₁₊_ir, L₂_ir, x_ir = diffusionoperators(irregularGrid, (Mixed(ξ = ξ_lb), Mixed(ξ = ξ_ub)))
        @test L₁₋bc ≈ L₁₋bc_ir
        @test L₁₊bc ≈ L₁₊_ir
        @test L₂bc ≈ L₂_ir
        @test x ≈ x_ir
    end
end

@testset "Operators under homogenous absorbing boundary conditions" begin
    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋bc = L₁₋bc(x̄, bc), L₁₊bc = L₁₊bc(x̄, bc), L₂bc = L₂bc(x̄, bc), x = interiornodes(x̄, bc))
    diffusionoperators_without_bc(x̄) = (L₁₋ = L₁₋(x̄), L₁₊ = L₁₊(x̄), L₂ = L₂(x̄), x = interiornodes(x̄))


    @testset "Accuracy & regression test" begin
        # Uniform grid
        x̄ = 0:4
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Absorbing(), Absorbing()))
        L₁₋, L₁₊, L₂, x = diffusionoperators_without_bc(x̄)
        
        # under absorbing barrier conditions, the corresponding 
        # operators on the interiors are identical as interior of extended operators
        @test L₁₋bc == L₁₋[:,2:end-1]
        @test L₁₊bc == L₁₊[:,2:end-1]
        @test L₂bc == L₂[:,2:end-1]

        # check accuracy as regression tests
        @test Array(L₁₋bc) == [1. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -1.]
        @test Array(L₂bc) == [-2. 1. 0.; 1. -2. 1.; 0. 1. -2.]
        @test Array(x) == [1; 2; 3]

        # Reflecting on lb, absorbing on ub
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Reflecting(), Absorbing()))
        L₁₋, L₁₊, L₂, x = diffusionoperators_without_bc(x̄)

        # check accuracy as regression tests
        @test Array(L₁₋bc) == [0. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1. 1.; 0. 0. -1.]
        @test Array(L₂bc) == [-1. 1. 0.; 1. -2. 1.; 0. 1. -2.]
        @test Array(x) == [1; 2; 3]
        
        # Absorbing on lb, reflecting on ub
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Absorbing(), Reflecting()))
        L₁₋, L₁₊, L₂, x = diffusionoperators_without_bc(x̄)

        # check accuracy as regression tests
        @test Array(L₁₋bc) == [1. 0. 0.; -1. 1. 0.; 0. -1. 1.]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1. 1.; 0. 0. 0.]
        @test Array(L₂bc) == [-2. 1. 0.; 1. -2. 1.; 0. 1. -1.]
        @test Array(x) == [1; 2; 3]

        # Irregular grid (with Δ_1m = Δ_1p and Δ_Mp = Δ_Mm)
        x̄ = [0.0; 1.0; 2.0; 4.0; 6.0]
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Absorbing(), Absorbing()))
        L₁₋, L₁₊, L₂, x = diffusionoperators_without_bc(x̄)
        
        # under absorbing barrier conditions, the corresponding 
        # operators on the interiors are identical as interior of extended operators
        @test L₁₋bc == L₁₋[:,2:end-1]
        @test L₁₊bc == L₁₊[:,2:end-1]
        @test L₂bc == L₂[:,2:end-1]

        # check accuracy as regression tests
        @test Array(L₁₋bc) == [1. 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. -1/2]
        @test Array(L₂bc) == [-2. 1. 0.; 2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 0. 1/4 -1/2]
        @test Array(x) == [1.0; 2.0; 4.0]

        # Irregular grid (with Δ_1m != Δ_1p and Δ_Mp != Δ_Mm)
        x̄ = [-1.0; 1.0; 2.0; 4.0; 7.0]

        Δ_1m = x̄[2] - x̄[1]
        Δ_1p = x̄[3] - x̄[2]
        Δ_Mm = x̄[end-1] - x̄[end-2]
        Δ_Mp = x̄[end] - x̄[end-1]
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(x̄, (Absorbing(), Absorbing()))
        L₁₋, L₁₊, L₂, x = diffusionoperators_without_bc(x̄)
        
        # under absorbing barrier conditions, the corresponding 
        # operators on the interiors are identical as interior of extended operators
        @test L₁₋bc == L₁₋[:,2:end-1]
        @test L₁₊bc == L₁₊[:,2:end-1]
        @test L₂bc == L₂[:,2:end-1]

        # check accuracy as regression tests
        @test Array(L₁₋bc) == [1/Δ_1m 0. 0.; -1. 1. 0.; 0. -1/2 1/2]
        @test Array(L₁₊bc) == [-1. 1. 0.; 0. -1/2 1/2; 0. 0. -1/Δ_Mp]
        Ξ_1 = -2*(1/(Δ_1m*Δ_1p)-1/((Δ_1p+Δ_1m)*(Δ_1m)))
        Ξ_M = -2*(1/(Δ_Mm*Δ_Mp)-1/((Δ_Mp+Δ_Mm)*(Δ_Mp)))
        
        @test Array(L₂bc) == [-2/(Δ_1p*Δ_1m) 2/(Δ_1p*(Δ_1m+Δ_1p)) 0.; 
                            2/((2.0+1.0)*1.0) -2/(2.0*1.0) 2/((2.0+1.0)*2.0); 
                            0. 2/(Δ_Mm*(Δ_Mm+Δ_Mp)) -2/(Δ_Mm*Δ_Mp) ]
        @test Array(x) == [1.0; 2.0; 4.0]
    end

    @testset "Consistency" begin
        # Test for consistency
        uniformGrid = range(0.0, 1.0, length = 500)
        irregularGrid = collect(uniformGrid)
        L₁₋bc, L₁₊bc, L₂bc, x = diffusionoperators(uniformGrid,(Absorbing(), Absorbing()))
        L₁₋bc_ir, L₁₊_ir, L₂_ir, x_ir = diffusionoperators(irregularGrid, (Absorbing(), Absorbing()))
        @test L₁₋bc ≈ L₁₋bc_ir
        @test L₁₊bc ≈ L₁₊_ir
        @test L₂bc ≈ L₂_ir
    end
end

@testset "Operators under upwind schemes" begin
    f(x) = x^2
    μ(x) = -x # drift depends on state
    σ = 0.1
    ρ = 0.05
    M = 100 # size of grid
    x̄ = range(-1.0, 1.0, length = M+2) # grid
    x = interiornodes(x̄)
    bc = (Reflecting(), Reflecting()) # specify BC (reflecting barrier)
    ## M-vector of drifts stacked according to the states
    μs = μ.(x)

    # Define first order differential operator using upwind scheme
    L_1_upwind = (μs .<= 0) .* L₁₋bc(x̄, bc) + (μs .> 0) .* L₁₊bc(x̄, bc)

    indices_m = findall(μs .<= 0)
    indices_p = findall(μs .> 0)

    @test L_1_upwind[indices_m, indices_m] == Array(L₁₋bc(x̄, bc)[indices_m, indices_m])
    @test L_1_upwind[indices_p, indices_p] == Array(L₁₊bc(x̄, bc)[indices_p, indices_p])
end

@testset "Input Type Variance" begin
    # helper function to get differential operators in a more efficient way
    diffusionoperators(x̄, bc) = (L₁₋bc = L₁₋bc(x̄, bc), L₁₊bc = L₁₊bc(x̄, bc), L₂bc = L₂bc(x̄, bc))

    # BigFloat
    uniformGrid = range(BigFloat(0.0), BigFloat(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Mixed(ξ = one(BigFloat)), Mixed(ξ = one(BigFloat))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(ξ = one(BigFloat)), Mixed(ξ = one(BigFloat))))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Mixed(ξ = one(BigFloat)), Mixed(ξ = one(BigFloat))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(ξ = one(BigFloat)), Mixed(ξ = one(BigFloat))))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)

    # Float32
    uniformGrid = range(Float32(0.0), Float32(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Mixed(ξ = one(Float32)), Mixed(ξ = one(Float32))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(ξ = one(Float32)), Mixed(ξ = one(Float32))))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc = L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Mixed(ξ = one(Float32)), Mixed(ξ = one(Float32))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(ξ = one(Float32)), Mixed(ξ = one(Float32))))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc = L₁₊_cache, L₂bc = L₂_cache)

    # Duals
    uniformGrid = range(Dual(0.0), Dual(1.0), length = 100)
    irregularGrid = collect(uniformGrid)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(uniformGrid, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))
    @test @inferred(diffusionoperators(irregularGrid, (Reflecting(), Reflecting()))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(uniformGrid, (Mixed(ξ = Dual(1.0)), Mixed(ξ = Dual(1.0))))
    @test @inferred(diffusionoperators(uniformGrid, (Mixed(ξ = Dual(1.0)), Mixed(ξ = Dual(1.0))))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
    L₁₋bc_cache, L₁₊_cache, L₂_cache = diffusionoperators(irregularGrid, (Mixed(ξ = Dual(1.0)), Mixed(ξ = Dual(1.0))))
    @test @inferred(diffusionoperators(irregularGrid, (Mixed(ξ = Dual(1.0)), Mixed(ξ = Dual(1.0))))) == (L₁₋bc = L₁₋bc_cache, L₁₊bc =  L₁₊_cache, L₂bc = L₂_cache)
end
