Lₙbc_convenience(x̄, bc) = Lₙbc(x̄, bc, JumpProcess(x̄, -1)) 
@testset "Accuracy test for 1 by 1 with two states" begin
    for L1_bc_generator in (L₁₊bc, L₁₋bc, L₂bc, Lₙbc_convenience),
        L2_bc_generator in (L₁₊bc, L₁₋bc, L₂bc, Lₙbc_convenience)
        M = 1 
        x̄ = 1:(M+2)
        bc = (Reflecting(), Reflecting())
        L1_bc = L1_bc_generator(x̄, bc)
        L2_bc = L2_bc_generator(x̄, bc)

        Q2 = zeros(2,2)

        @test [L1_bc[1,1] 0; 0 L2_bc[1,1]] == jointoperator_bc((L1_bc, L2_bc), Q2)
    end
end
@testset "Accuracy test for 2 by 2 with two states" begin
    for L1_bc_generator in (L₁₊bc, L₁₋bc, L₂bc, Lₙbc_convenience),
        L2_bc_generator in (L₁₊bc, L₁₋bc, L₂bc, Lₙbc_convenience),

        M = 2
        x̄ = 1:(M+2)
        bc = (Reflecting(), Reflecting())
        L1_bc = L1_bc_generator(x̄, bc)
        L2_bc = L2_bc_generator(x̄, bc)

        Q2 = zeros(2,2)

        @test [L1_bc[1,1] L1_bc[1,2] 0 0; 
               L1_bc[2,1] L1_bc[2,2] 0 0;
               0 0 L2_bc[1,1] L2_bc[1,2];
               0 0 L2_bc[2,1] L2_bc[2,2]] == jointoperator_bc((L1_bc, L2_bc), Q2)
    end
end