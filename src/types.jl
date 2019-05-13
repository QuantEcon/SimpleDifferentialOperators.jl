# Types for Boundary Conditions
    # Abstracts
    abstract type BoundaryCondition end
    abstract type HomogeneousBoundaryCondition <: BoundaryCondition end

# Types for Differences
    abstract type DiscretizationMethod{N} end
    struct DefaultFirstDifference <: DiscretizationMethod{1} end
    struct ForwardFirstDifference <: DiscretizationMethod{1} end
    struct BackwardFirstDifference <: DiscretizationMethod{1} end
    struct CentralSecondDifference <: DiscretizationMethod{2} end
    struct JumpProcess{T} where T <: AbstractArray
        jumps::T
    end

# Concretes
    struct Reflecting <: HomogeneousBoundaryCondition end
    @with_kw struct Mixed{T} <: HomogeneousBoundaryCondition where T <: Real
        ξ::T = 0.0
        direction = :auto # :forward/:backward/:auto
    end
    struct Absorbing <: HomogeneousBoundaryCondition end