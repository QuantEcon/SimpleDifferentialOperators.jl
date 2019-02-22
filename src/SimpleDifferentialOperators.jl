module SimpleDifferentialOperators

# Pkg Dependencies
using LinearAlgebra, SparseArrays

# Includes
include("types.jl")
include("operators.jl")
include("operators-without-bc.jl")
include("basis.jl")

# Exports
# Types
export BoundaryCondition,
       HomogeneousBoundaryCondition,
       InhomogeneousBoundaryCondition,
       Reflecting,
       Mixed,
       Absorbing,
       NoBoundary

# Functions
export diffusionoperators

end # module
