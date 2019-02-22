using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@testset "Operators without boundary conditions" begin
    # Check consistency
    for N in 3:10 # repeat with different node numbers
        # setup 
        f(x) = x^2
        μ = -0.1 # constant negative drift
        σ = 0.1
        N = 3
        x = range(0.0, 1.0, length = N)

        # operators with reflecting boundary conditions
        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())
        A = μ*L_1_minus + σ^2 / 2 * L_2 
        ## solve the value function
        v_bc = (I * 0.05 - A) \ f.(x) 

        # operators without boundary conditions, adding extra two rows for boundary conditions
        L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x)
        L = μ*L_1_minus + σ^2 / 2 * L_2 
        B = transpose([[-1; 1; zeros(N)] [zeros(N); -1; 1]])
        A = [([zeros(N) Diagonal(ones(N,N)) zeros(N)] * 0.05 - L); B]

        ## solve the value function including the boundary conditions
        v_bar = A \ [f.(x); 0.0; 0.0] 
        v_interior = v_bar[2:end-1] # extract the interior

        @test v_bc ≈ v_interior
    end
end