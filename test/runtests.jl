using SimpleDifferentialOperators
using Test, LinearAlgebra, PATHSolver, Suppressor, Parameters, Random, BandedMatrices
using DualNumbers

@elapsed begin
    @time @testset "Operators with boundary conditions" begin include("operators-fd.jl") end
    @time @testset "Operators without boundary conditions" begin include("operators-fd-without-bc.jl") end
    @time @testset "Operators without boundary conditions for jump diffusion" begin include("operators-jump.jl") end 
    @time @testset "Linear Complementarity Problems" begin include("lcp.jl") end
    @time @testset "Boundary extrapolation" begin include("utilities/extrapolatetoboundary.jl") end
    @time @testset "Internal helper functions" begin include("utilities/findnearestindex.jl") end
end
