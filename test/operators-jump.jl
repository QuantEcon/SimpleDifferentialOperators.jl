using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers


@testset "Truncation rules when indices are explicitly given" begin
    x̄ = 1:5
    jumps = [0; -1; 0]
    @test JumpProcess(x̄, jumps).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps 

    jumps = [0; 1; 0]
    @test JumpProcess(x̄, jumps).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps 

    jumps = [0; -2; 0]
    @test JumpProcess(x̄, jumps).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps 
    
    jumps = [0; 2; 0]
    @test JumpProcess(x̄, jumps).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == jumps
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == [0; 1; 0] 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps
    
    jumps = [-1; -2; 0]
    @test JumpProcess(x̄, jumps).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps 
    
    jumps = [0; 2; 1]
    @test JumpProcess(x̄, jumps).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == jumps
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == [0; 1; 0] 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps
    
    jumps = [-1; -2; 1]
    @test JumpProcess(x̄, jumps).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == [0; -1; 1]
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == [-1; -2; 0]
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps
    
    jumps = [-1; 2; 1]
    @test JumpProcess(x̄, jumps).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumps, (:interior, :boundary)).jumps == [0; 2; 1]
    @test JumpProcess(x̄, jumps, (:boundary, :interior)).jumps == [-1; 1; 0] 
    @test JumpProcess(x̄, jumps, (:boundary, :boundary)).jumps == jumps
end

@testset "Truncation rules when jump functions are provided" begin
    x̄ = 1:5
    jumps = [0; -1; 0]
    jumpf = x -> (x == 3) ? -1 : 0 
    @test JumpProcess(x̄, jumpf).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps 

    jumps = [0; 1; 0]
    jumpf = x -> (x == 3) ? 1 : 0 
    @test JumpProcess(x̄, jumpf).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps 

    jumps = [0; -2; 0]
    jumpf = x -> (x == 3) ? -2 : 0 
    @test JumpProcess(x̄, jumpf).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps 
    
    jumps = [0; 2; 0]
    jumpf = x -> (x == 3) ? 2 : 0 
    @test JumpProcess(x̄, jumpf).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == jumps
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == [0; 1; 0] 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps
    
    jumps = [-1; -2; 0]
    jumpf = x -> (x <= 3) ? -(x-1) : 0 
    @test JumpProcess(x̄, jumpf).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == jumps 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps 
    
    jumps = [0; 2; 1]
    jumpf = x -> (x >= 3) ? (5-x) : 0 
    @test JumpProcess(x̄, jumpf).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == jumps
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == [0; 1; 0] 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps
    
    jumps = [-1; -2; 1]
    jumpf = x -> (x <= 3) ? -(x-1) : 1
    @test JumpProcess(x̄, jumpf).jumps == [0; -1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == [0; -1; 1]
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == [-1; -2; 0]
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps
    
    jumps = [-1; 2; 1]
    jumpf = x -> (x >= 3) ? (5-x) : -1
    @test JumpProcess(x̄, jumpf).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == [0; 2; 1]
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == [-1; 1; 0] 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps

    # on irregular grids
    x̄ = [1; 2.5; 3; 4.9; 5]
    jumps = [-1; 2; 1]
    jumpf = x -> (x >= 3) ? (5-x) : -1
    @test JumpProcess(x̄, jumpf).jumps == [0; 1; 0]
    @test JumpProcess(x̄, jumpf, (:interior, :boundary)).jumps == [0; 2; 1]
    @test JumpProcess(x̄, jumpf, (:boundary, :interior)).jumps == [-1; 1; 0] 
    @test JumpProcess(x̄, jumpf, (:boundary, :boundary)).jumps == jumps
end