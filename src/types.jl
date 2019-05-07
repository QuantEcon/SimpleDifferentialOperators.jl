# Types for Boundary Conditions
    # Abstracts
    abstract type BoundaryCondition end
    abstract type HomogeneousBoundaryCondition <: BoundaryCondition end
    abstract type InhomogeneousBoundaryCondition <: BoundaryCondition end

# Types for Differences
    abstract type DifferenceMethod{N} end
    struct DefaultFirstDifference <: DifferenceMethod{1} end
    struct ForwardFirstDifference <: DifferenceMethod{1} end
    struct BackwardFirstDifference <: DifferenceMethod{1} end
    struct CentralSecondDifference <: DifferenceMethod{2} end

# Concretes
    struct Reflecting <: HomogeneousBoundaryCondition end
    @with_kw struct Mixed{T} <: HomogeneousBoundaryCondition where {T <: DifferenceMethod{1}}
        ξ::Real = 0.0
        direction::T = DefaultFirstDifference()
    end