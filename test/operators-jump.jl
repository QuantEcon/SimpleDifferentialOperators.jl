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

@testset "Accuracy test for jump diffusion extension operators" begin
    x̄ = [1; 2.3; 3.4; 4.5; 5.0]
    jumps = [0; -1; 0]
    method = JumpProcess(x̄, jumps)
    L̄ = ExtensionDifferentialOperator(x̄, method) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. 0. 0.]

    jumps = [0; 1; 0]
    method = JumpProcess(x̄, jumps)
    L̄ = ExtensionDifferentialOperator(x̄, method) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]

    jumps = [0; -2; 0]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 1. 0. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 1. 0. -1. 0. 0.; 0. 0. 0. 0. 0.]

    jumps = [0; 2; 0]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 0. 1.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 0. 1.; 0. 0. 0. 0. 0.]

    jumps = [-1; -2; 0]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 1. 0. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 1. 0. -1. 0. 0.; 0. 0. 0. 0. 0.]

    jumps = [0; 2; 1]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 0. 1.; 0. 0. 0. -1. 1.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 0. 1.; 0. 0. 0. -1. 1.]
    
    jumps = [-1; -2; 1]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 0. -1. 1.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 1. 0. -1. 0. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 1. 0. -1. 0. 0.; 0. 0. 0. -1. 1.]

    jumps = [-1; 2; 1]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 0. -1. 0. 1.; 0. 0. 0. -1. 1.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 0. 0. -1. 0. 1.; 0. 0. 0. -1. 1.]

    jumps = [-1; -1; -1]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 1. -1. 0.]
    @test L̄ == ExtensionDifferentialOperator(x̄, JumpProcess(x̄, -1))
    @test L̄ == ExtensionDifferentialOperator(x̄, JumpProcess(x̄, -1.0))
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. 0. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 1. -1. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 1. -1. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [1. -1. 0. 0. 0.; 0. 1. -1. 0. 0.; 0. 0. 1. -1. 0.]

    jumps = [1; 1; 1]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps)) 
    @test L̄ == ExtensionDifferentialOperator(x̄, JumpProcess(x̄, 1))
    @test L̄ == ExtensionDifferentialOperator(x̄, JumpProcess(x̄, 1.0))
    @test Array(L̄) == [0. -1. 1. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L̄) == [0. -1. 1. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. -1. 1.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L̄) == [0. -1. 1. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. 0. 0.]
    L̄ = ExtensionDifferentialOperator(x̄, JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L̄) == [0. -1. 1. 0. 0.; 0. 0. -1. 1. 0.; 0. 0. 0. -1. 1.]
end


@testset "Accuracy test for jump diffusion operators on the interior" begin
    x̄ = [1; 2.3; 3.4; 4.5; 5.0]
    jumps = [0; -1; 0]
    bc = (Absorbing(), Absorbing())
    method = JumpProcess(x̄, jumps)
    L = Lₙbc(x̄,bc,method) 
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0.]

    jumps = [0; 1; 0]
    method = JumpProcess(x̄, jumps)
    L = Lₙbc(x̄,bc,method) 
    @test Array(L) == [0. 0. 0.; 0. -1. 1.; 0. 0. 0.]

    jumps = [0; -1; 0]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps))
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]

    jumps = [0; -2; 0]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps))
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:boundary, :interior))) 
    @test Array(L) == [0. 0. 0.; 0. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:boundary, :boundary))) 
    @test Array(L) == [0. 0. 0.; 0. -1. 0.; 0. 0. 0. ]

    bc = (Reflecting(), Reflecting())
    jumps = [0; -2; 0]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps))
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    L = Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:interior, :boundary))) 
    @test Array(L) == [0. 0. 0.; 1. -1. 0.; 0. 0. 0. ]
    @test_throws ErrorException Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:boundary, :interior)))  
    @test_throws ErrorException Lₙbc(x̄,bc,JumpProcess(x̄, jumps, (:boundary, :boundary))) 
end