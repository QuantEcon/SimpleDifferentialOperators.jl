using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Accuracy & regression test" begin
    regular_grid = 0:3
    irregular_grid = [0.0; 1.0; 2.0; 4.0; 6.0]
    ## regular grids
    L₁₋_regular = L₁₋(regular_grid)
    L̄₁₊_regular = L₁₊(regular_grid)
    L̄₂_regular = L₂(regular_grid)
    x_regular = interiornodes(regular_grid)
    @test Array(L₁₋_regular) == [-1. 1. 0. 0.; 0. -1. 1. 0.]
    @test Array(L̄₁₊_regular) == [0. -1. 1. 0.; 0. 0. -1. 1.]
    @test Array(L̄₂_regular) == [1. -2. 1. 0.; 0. 1. -2. 1.]
    @test Array(x_regular) == [1; 2]

    @test @inferred(L₁₋(regular_grid)) == L₁₋_regular
    @test @inferred(L₁₊(regular_grid)) == L̄₁₊_regular
    @test @inferred(L₂(regular_grid)) == L̄₂_regular
    @test @inferred(interiornodes(regular_grid)) == x_regular

    ## irregular grids
    L₁₋_irregular = L₁₋(irregular_grid)
    L̄₁₊_irregular = L₁₊(irregular_grid)
    L̄₂_irregular = L₂(irregular_grid)
    x_irregular = interiornodes(irregular_grid)
    @test Array(L₁₋_irregular) == [-1. 1. 0. 0. 0.; 0. -1. 1. 0. 0.; 0. 0. -1/2 1/2 0.]
    @test Array(L̄₁₊_irregular) == [0. -1. 1. 0. 0.; 0. 0. -1/2 1/2 0.; 0. 0. 0. -1/2 1/2]
    @test Array(L̄₂_irregular) == [1. -2. 1. 0. 0.; 0 2/((1+2)*1) -2/(2*1) 2/((1+2)*2) 0; 0. 0. 1/4 -2/(2*2) 1/4]
    @test Array(x_irregular) == [1; 2; 4]

    @test @inferred(L₁₋(irregular_grid)) == L₁₋_irregular
    @test @inferred(L₁₊(irregular_grid)) == L̄₁₊_irregular
    @test @inferred(L₂(irregular_grid)) == L̄₂_irregular
    @test @inferred(interiornodes(irregular_grid)) == x_irregular
end

@testset "Consistency" begin
    regular_grids = [range(0.0, 1.0, length = 5), range(-1.0, 1.0, length = 5)]
    irregular_grids = [[-1.0; 0.0; 2.0; 5.0; 9.0], [0.1; 0.3; 0.9; 1.9; 2.5; 9.0]]
    # Check consistency
    for x̄ in [regular_grids; irregular_grids]
        # setup
        f(x) = x^2
        μ = -0.1 # constant negative drift
        σ = 0.1
        ρ = 0.05
        
        x = interiornodes(x̄)
        M = length(x)
        bc = (Reflecting(), Reflecting()) # specify BC (reflecting barrier)

        # operators with reflecting boundary conditions
        L = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)
        ## solve the value function
        v = (I * ρ - L) \ f.(x)

        # operators without boundary conditions, adding extra two rows for boundary conditions
        L̄ = μ*L₁₋(x̄) + σ^2 / 2 * L₂(x̄)
        B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])
        L = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * ρ - L̄); B]

        ## solve the value function including the boundary conditions
        v̄ = L \ [f.(x); 0.0; 0.0]

        # currently broken on regular grids
        @test v ≈ v̄[2:end-1]
    end
end
