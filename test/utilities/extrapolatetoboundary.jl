
@testset "Operators under mixed barrier conditions" begin
    @testset "Accuracy & regression test" begin
    regular_grids = [0:4, -1:6]
    irregular_grids = [[-1.0; 0.0; 2.0; 5.0; 9.0], [0.2; 0.3; 0.5; 0.9; 1.6; 2.5]]
        for x̄ in [regular_grids; irregular_grids], 
            (L_bc_generator, L_generator) in [(L₁₋bc, L₁₋), (L₁₊bc, L₁₊), (L₂bc, L₂)],
            (ξ_lb, ξ_ub) in [(-2, 1), (-3, -2), (3, -4), (4, 3)]
            # Setup
            # RHS function
            f(x) = x^2

            # grid setup
            # x̄ = 0:4
            x = interiornodes(x̄)
            M = length(x)
            Δ_0p = x̄[2] - x̄[1] 
            Δ_Mp = x̄[end] - x̄[end-1]

            # boundary condition
            bc = (Mixed(ξ_lb), Mixed(ξ_ub))

            # corresponding boundary condition matrix
            B = transpose([[-1+ξ_lb*Δ_0p; 1; zeros(M)] [zeros(M); -1; 1+ξ_ub*Δ_Mp]])
            b = [0; 0]

            # operator on the interior, with bc applied 
            L_bc = L_bc_generator(x̄, bc)
            # operator on the exterior
            L = L_generator(x̄)

            # solution for interior v from L_bc
            v_bc = L_bc \ f.(x)
            # solution for extended v from L by stacking up boundary condition matrices
            v̄ = [L; B] \ [f.(x); b]

            # test if extrapolated solutions are identical
            @test extrapolatetoboundary(v_bc, x̄, bc) ≈ v̄ 
            
        end
    end
end