using SimpleDifferentialOperators
using Test, LinearAlgebra

@testset "SimpleDifferentialOperators.jl" begin
    grid = range(1., 10., length = 100) # calls the regular method
    @test which(rescaled_diffusionoperators, (typeof(grid), typeof(0.))).sig == Tuple{typeof(rescaled_diffusionoperators),AbstractRange,Any}
    grid = collect(grid) # will call the irregular method
    @test which(rescaled_diffusionoperators,  (typeof(grid), typeof(0.))).sig == Tuple{typeof(rescaled_diffusionoperators),AbstractArray,Any}
end
