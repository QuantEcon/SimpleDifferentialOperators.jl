using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers, PATHSolver, Suppressor, Parameters

@elapsed begin
    @time @testset "Operators with boundary conditions" begin include("operators.jl") end
    @time @testset "Operators without boundary conditions" begin include("operators-without-bc.jl") end
    @time @testset "Linear Complementarity Problems" begin include("lcp.jl") end
end
