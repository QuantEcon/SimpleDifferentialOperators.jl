using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Accuracy & regression test" begin
    uniform_grid = 1:1:2
    irregular_grid = collect(uniform_grid)
    ## regular grids
    L̄₁₋, L̄₁₊, L̄₂, x̄ = diffusionoperators(uniform_grid)
    @test Array(L̄₁₋) == [-1. 1. 0. 0.; 0. -1. 1. 0.]
    @test Array(L̄₁₊) == [0. -1. 1. 0.; 0. 0. -1. 1.]
    @test Array(L̄₂) == [1. -2. 1. 0.; 0. 1. -2. 1.]
    @test Array(x̄) == [0; 1; 2; 3]
    @test @inferred(diffusionoperators(uniform_grid)) == (L̄₁₋ = L̄₁₋, L̄₁₊ = L̄₁₊, L̄₂ = L̄₂, x̄ = x̄)

    ## irregular grids
    L̄₁₋, L̄₁₊, L̄₂, x̄ = diffusionoperators(irregular_grid)
    @test Array(L̄₁₋) == [-1. 1. 0. 0.; 0. -1. 1. 0.]
    @test Array(L̄₁₊) == [0. -1. 1. 0.; 0. 0. -1. 1.]
    @test Array(L̄₂) == [1. -2. 1. 0.; 0. 1. -2. 1.]
    @test Array(x̄) == [0; 1; 2; 3]
    @test @inferred(diffusionoperators(irregular_grid)) == (L̄₁₋ = L̄₁₋, L̄₁₊ = L̄₁₊, L̄₂ = L̄₂, x̄ = x̄)
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
        x = range(0.0, 1.0, length = N)

        # operators with reflecting boundary conditions
        L₁₋, L₁₊, L₂, x̄ = diffusionoperators(x, (Reflecting(), Reflecting()))
        A = μ*L₁₋ + σ^2 / 2 * L₂
        ## solve the value function
        v_bc = (I * ρ - A) \ f.(x)

        # operators without boundary conditions, adding extra two rows for boundary conditions
        L̄₁₋, L̄₁₊, L̄₂, x̄ = diffusionoperators(x)
        L = μ*L̄₁₋ + σ^2 / 2 * L̄₂
        B = transpose([[-1; 1; zeros(N)] [zeros(N); -1; 1]])
        A = [([zeros(N) Diagonal(ones(N,N)) zeros(N)] * 0.05 - L); B]

        ## solve the value function including the boundary conditions
        v_bar = A \ [f.(x); 0.0; 0.0]
        v_interior = v_bar[2:end-1] # extract the interior

        @test v_bc ≈ v_interior
    end
end
