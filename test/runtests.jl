using SimpleDifferentialOperators
using Test, LinearAlgebra

@testset "Dispatch Tests" begin
    grids = [range(1.0, 10.0, length = 100), collect(range(1.0, 10.0, length = 100))]
    BCs = [Reflecting(), Mixed(2.0), Absorbing(1.0, 2.0)] # list of all possible BCs
    inputMatrix = collect(Base.product(grids, BCs, BCs))
    for input in inputMatrix
        res = diffusionoperators(input...)
        # Test grid dispatch
        if input[1] isa AbstractRange
            @test res[1] == "Uniform"
        elseif input[1] isa AbstractArray
            @test res[1] == "Irregular"
        else
            @error "This should never happen."
        end

        # Test BC1
        if input[2] isa Reflecting
            @test res[2] == "Reflecting"
        elseif input[2] isa Mixed
            @test res[2] == "Mixed"
        elseif input[2] isa Absorbing
            @test res[2] == "Absorbing"
        else
            @error "This should never happen."
        end

        # Test BC2
        if input[3] isa Reflecting
            @test res[3] == "Reflecting"
        elseif input[3] isa Mixed
            @test res[3] == "Mixed"
        elseif input[3] isa Absorbing
            @test res[3] == "Absorbing"
        else
            @error "This should never happen."
        end
    end
end
