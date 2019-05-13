using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers


@testset "Truncation rules" begin
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