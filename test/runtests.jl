using SimpleDifferentialOperators
using Test, LinearAlgebra, PATHSolver, Suppressor, Parameters, Random, BandedMatrices
using DualNumbers

@elapsed begin
    @time @testset "Jump diffusion" begin include("jump-diffusion.jl") end 
    @time @testset "Operators with boundary conditions" begin include("operators.jl") end
    @time @testset "Operators without boundary conditions" begin include("operators-without-bc.jl") end
    @time @testset "Linear Complementarity Problems" begin include("lcp.jl") end
    @time @testset "Boundary extrapolation" begin include("utilities/extrapolatetoboundary.jl") end
    @time @testset "Internal helper functions" begin include("utilities/findnearestindex.jl") end
end
