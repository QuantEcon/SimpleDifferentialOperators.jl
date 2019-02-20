module SimpleDifferentialOperators

# Pkg Dependencies
using LinearAlgebra

# Includes
include("types.jl")
include("operators.jl")
include("basis.jl")

# Exports
# Types
export BoundaryCondition,
       HomogeneousBoundaryCondition,
       InhomogeneousBoundaryCondition,
       Reflecting,
       Mixed,
       Absorbing

# Functions
export diffusionoperators

end # module
