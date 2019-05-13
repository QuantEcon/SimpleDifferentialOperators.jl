# Types for Boundary Conditions
    # Abstracts
    abstract type BoundaryCondition end
    abstract type HomogeneousBoundaryCondition <: BoundaryCondition end

# Types for Differences
    abstract type DiscretizationMethod end
    struct DefaultFirstDifference <: DiscretizationMethod end
    struct ForwardFirstDifference <: DiscretizationMethod end
    struct BackwardFirstDifference <: DiscretizationMethod end
    struct CentralSecondDifference <: DiscretizationMethod end
    struct JumpProcess{T} <: DiscretizationMethod where T <: AbstractArray
        jumps::T
    end

# Concretes
    struct Reflecting <: HomogeneousBoundaryCondition end
    @with_kw struct Mixed{T} <: HomogeneousBoundaryCondition where T <: Real
        ξ::T = 0.0
        direction = :auto # :forward/:backward/:auto
    end
    struct Absorbing <: HomogeneousBoundaryCondition end