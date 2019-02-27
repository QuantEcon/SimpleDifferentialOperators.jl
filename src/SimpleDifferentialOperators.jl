module SimpleDifferentialOperators

# Pkg Dependencies
using LinearAlgebra, SparseArrays

# Includes
include("types.jl")
include("basis.jl")
include("operators.jl")

# Exports
# Boundary Conditions
export BoundaryCondition,
       HomogeneousBoundaryCondition,
       InhomogeneousBoundaryCondition,
       Reflecting,
       Mixed,
       Absorbing,
       NoBoundary

# Differential Types
export DifferenceMethod,
        ForwardFirstDifference,
        BackwardFirstDifference,
        CentralSecondDifference

# Functions
export DifferentialOperator,
        L₁₋,
        L₁₊,
        L₂,
        x̄,
        diffusionoperators,
        L̄₁₊,
        L̄₁₋,
        L̄₂

end # module
