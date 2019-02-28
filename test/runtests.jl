using SimpleDifferentialOperators
using Test, LinearAlgebra, DualNumbers

@elapsed begin
    @time @testset "Operators with boundary conditions" begin include("operators.jl") end
    @time @testset "Operators without boundary conditions" begin include("operators-without-bc.jl") end
end
