module SimpleDifferentialOperators

# Pkg Dependencies
using LinearAlgebra

# Includes
include("operators.jl")

# Exports
export robin_diffusionoperators, # both the irregular and regular method
       reflecting_diffusionoperators # with ξ = 0

end # module
