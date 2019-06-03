
@testset "Operators under mixed barrier conditions" begin
    regular_grids = [0:4, -1:6]
    irregular_grids = [[-1.0; 0.0; 2.0; 5.0; 9.0], [0.2; 0.3; 0.5; 0.9; 1.6; 2.5]]
        for x̄ in [regular_grids; irregular_grids], 
            (L_bc_generator, L_generator) in [(L₁₋bc, L₁₋), (L₁₊bc, L₁₊), (L₂bc, L₂)],
            (ξ_lb, ξ_ub) in [(-2, 1), (-3, -2), (3, -4), (4, 3)],
            directions in [(:forward, :backward), (:backward, :forward),
                            (:backward, :backward), (:forward, :forward)]
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
            bc = (Mixed(ξ = ξ_lb, direction = directions[1]), 
                Mixed(ξ = ξ_ub, direction = directions[2]))

            # corresponding boundary condition matrix
            B_lb = (directions[1] != :backward) ? [-1+ξ_lb*Δ_0p; 1; zeros(M)] : [-1; 1+ξ_lb*Δ_0p; zeros(M)]
            B_ub = (directions[2] != :forward) ? [zeros(M); -1; 1+ξ_ub*Δ_Mp] : [zeros(M); -1+ξ_ub*Δ_Mp; 1]

            B = transpose([B_lb B_ub])
            b = [0; 0]

            # operator on the interior, with bc applied 
            L_bc = L_bc_generator(x̄, bc)
            # operator on the exterior
            L = L_generator(x̄)

            # solution for interior v from L_bc
            v_bc = L_bc \ f.(x)
            # solution for extended v from L by stacking up boundary condition matrices
            v̄ = [L; B] \ [f.(x); b]

            # test if interior solutions are identical
            @test v_bc ≈ v̄[2:end-1]
            # test if extrapolated solutions are identical
            @test extrapolatetoboundary(v_bc, x̄, bc) ≈ v̄ 
        end
end

@testset "Accuracy test for absorbing boundary conditions" begin
    regular_grids = [0:4, -1:6]
    irregular_grids = [[-1.0; 0.0; 2.0; 5.0; 9.0], [0.2; 0.3; 0.5; 0.9; 1.6; 2.5]]
    for x̄ in [regular_grids; irregular_grids], 
        (L_bc_generator, L_generator) in [(L₁₋bc, L₁₋), (L₁₊bc, L₁₊), (L₂bc, L₂)]
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
        bc = (Absorbing(), Absorbing())

        # corresponding boundary condition matrix
        B = transpose([[1; zeros(M+1)] [zeros(M+1); 1]])
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

@testset "Accuracy test for nonhomogeneous absorbing boundary conditions" begin
    regular_grids = [0:4, -1:6]
    irregular_grids = [[-1.0; 0.0; 2.0; 5.0; 9.0], [0.2; 0.3; 0.5; 0.9; 1.6; 2.5]]
    for x̄ in [regular_grids; irregular_grids], 
        (L_bc_generator, L_generator, L_affine_generator) in [(L₁₋bc, L₁₋, L₁₋affine), (L₁₊bc, L₁₊, L₁₊affine), (L₂bc, L₂, L₂affine)],
        (S_lower, S_upper) in [(0.0, 1.0), (1.0, 0.0), (1.5, -1.5), (1.2, 1.3), (1.0, 1.0)]
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
        bc = (NonhomogeneousAbsorbing(S_lower), NonhomogeneousAbsorbing(S_upper))

        # corresponding boundary condition matrix
        B = transpose([[1; zeros(M+1)] [zeros(M+1); 1]])
        b = [S_lower; S_upper]

        # operator on the interior, with bc applied 
        L_bc = L_bc_generator(x̄, bc)
        # operator on the exterior
        L = L_generator(x̄)

        # solution for interior v from L_bc
        v_bc = L_bc \ (f.(x) + L_affine_generator(x̄, f.(x), bc))
        # solution for extended v from L by stacking up boundary condition matrices
        v̄ = [L; B] \ [f.(x); b]

        # test if extrapolated solutions are identical
        @test extrapolatetoboundary(v_bc, x̄, bc) ≈ v̄ 
        
    end
end