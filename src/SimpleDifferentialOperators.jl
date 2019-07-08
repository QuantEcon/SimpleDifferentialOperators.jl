module SimpleDifferentialOperators

# Pkg Dependencies
using LinearAlgebra, SparseArrays, Parameters, BandedMatrices, BlockBandedMatrices, LazyArrays

# Includes
include("types.jl")
include("operators.jl")
include("extensionoperators.jl")
include("jointoperators.jl")
include("affine.jl")
include("utilities/extrapolatetoboundary.jl")
include("utilities/findnearestindex.jl")

# Exports
# Boundary Conditions
export BoundaryCondition,
       HomogeneousBoundaryCondition,
       NonhomogeneousBoundaryCondition,
       Reflecting,
       Mixed,
       Absorbing,
       NonhomogeneousAbsorbing

# Differential Types
export DiscretizationMethod,
        ForwardFirstDifference,
        BackwardFirstDifference,
        CentralSecondDifference,
        JumpProcess

# Functions
export DifferentialOperator,
        ExtensionDifferentialOperator,
        L₁₋bc,
        L₁₊bc,
        L₂bc,
        Lₙbc,
        interiornodes,
        L₁₊,
        L₁₋,
        L₂,
        Lₙ,
        jointoperator_bc,
        extrapolatetoboundary,
        findnearestindex,
        L₀affine,
        L₁₋affine,
        L₁₊affine,
        L₂affine,
        Laffine

end # module
