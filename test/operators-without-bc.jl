using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Accuracy & regression test" begin
    regular_grid = 0:3
    irregular_grid = [0.0; 1.0; 2.0; 4.0; 6.0]
    ## regular grids
    L̄₁₋_regular = L̄₁₋(regular_grid)
    L̄₁₊_regular = L̄₁₊(regular_grid)
    L̄₂_regular = L̄₂(regular_grid)
    x_regular = interiornodes(regular_grid)
    @test Array(L̄₁₋_regular) == [-1. 1. 0. 0.; 0. -1. 1. 0.]
    @test Array(L̄₁₊_regular) == [0. -1. 1. 0.; 0. 0. -1. 1.]
    @test Array(L̄₂_regular) == [1. -2. 1. 0.; 0. 1. -2. 1.]
    @test Array(x_regular) == [1; 2]

    @test @inferred(L̄₁₋(regular_grid)) == L̄₁₋_regular
    @test @inferred(L̄₁₊(regular_grid)) == L̄₁₊_regular
    @test @inferred(L̄₂(regular_grid)) == L̄₂_regular
    @test @inferred(interiornodes(regular_grid)) == x_regular

    ## irregular grids
    L̄₁₋_irregular = L̄₁₋(irregular_grid)
    L̄₁₊_irregular = L̄₁₊(irregular_grid)
    L̄₂_irregular = L̄₂(irregular_grid)
    x_irregular = interiornodes(irregular_grid)
    @test Array(L̄₁₋_irregular) == [-1. 1. 0. 0. 0.; 0. -1. 1. 0. 0.; 0. 0. -1/2 1/2 0.]
    @test Array(L̄₁₊_irregular) == [0. -1. 1. 0. 0.; 0. 0. -1/2 1/2 0.; 0. 0. 0. -1/2 1/2]
    @test Array(L̄₂_irregular) == [1. -2. 1. 0. 0.; 0 2/((1+2)*1) -2/(2*1) 2/((1+2)*2) 0; 0. 0. 1/4 -2/(2*2) 1/4]
    @test Array(x_irregular) == [1; 2; 4]

    @test @inferred(L̄₁₋(irregular_grid)) == L̄₁₋_irregular
    @test @inferred(L̄₁₊(irregular_grid)) == L̄₁₊_irregular
    @test @inferred(L̄₂(irregular_grid)) == L̄₂_irregular
    @test @inferred(interiornodes(irregular_grid)) == x_irregular
end

@testset "Consistency" begin
    # Check consistency
    for N in 3:10 # repeat with different node numbers
        # setup
        f(x) = x^2
        μ = -0.1 # constant negative drift
        σ = 0.1
        ρ = 0.05
        N = 3
        x̄ = range(0.0, 1.0, length = N + 2)
        x = x̄[2:end-1]
        bc = (Reflecting(), Reflecting()) # specify BC (reflecting barrier)

        # operators with reflecting boundary conditions
        L = μ*L₁₋(x̄, bc) + σ^2 / 2 * L₂(x̄, bc)
        ## solve the value function
        v = (I * ρ - L) \ f.(x)

        # operators without boundary conditions, adding extra two rows for boundary conditions
        L̄ = μ*L̄₁₋(x̄) + σ^2 / 2 * L̄₂(x̄)
        B = transpose([[-1; 1; zeros(N)] [zeros(N); -1; 1]])
        L = [([zeros(N) Diagonal(ones(N,N)) zeros(N)] * ρ - L̄); B]

        ## solve the value function including the boundary conditions
        v̄ = L \ [f.(x); 0.0; 0.0]
        v_by_extended_operator = v̄[2:end-1] # extract the interior

        @test v ≈ v_by_extended_operator
    end
end
