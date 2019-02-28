# Types for Boundary Conditions
    # Abstracts
    abstract type BoundaryCondition end
    abstract type HomogeneousBoundaryCondition <: BoundaryCondition end
    abstract type InhomogeneousBoundaryCondition <: BoundaryCondition end
    # Concretes
    struct Reflecting <: HomogeneousBoundaryCondition end
    struct Mixed{T} <: HomogeneousBoundaryCondition where {T <: Real}
        ξ::T
    end

# Types for Differences
    abstract type DifferenceMethod{N} end
    struct ForwardFirstDifference <: DifferenceMethod{1} end
    struct BackwardFirstDifference <: DifferenceMethod{1} end
    struct CentralSecondDifference <: DifferenceMethod{2} end
