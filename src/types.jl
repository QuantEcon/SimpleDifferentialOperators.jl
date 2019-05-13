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
    function JumpProcess(x̄::AbstractArray, jumps::AbstractArray, 
                        truncate = (:interior, :interior))
        # each ith element of jumps defines the jump direction (negative/positive) and size
        # where i ranges from 1 to M (interior points)
        @assert length(x̄) == length(jumps) + 2 
        M = length(x̄) - 2

        # define default destinations and truncation indices
        destinations = Array(1:M) + jumps
        cutoff_lb = (truncate[1] == :interior) ? 1 : 0 
        cutoff_ub = (truncate[2] == :interior) ? M : M+1

        # perform truncations
        destinations[destinations .< cutoff_lb] .= cutoff_lb
        destinations[destinations .> cutoff_ub] .= cutoff_ub
        jumps_after_truncation = destinations - Array(1:M)

        return JumpProcess(jumps_after_truncation)
    end

# Concretes
    struct Reflecting <: HomogeneousBoundaryCondition end
    @with_kw struct Mixed{T} <: HomogeneousBoundaryCondition where T <: Real
        ξ::T = 0.0
        direction = :auto # :forward/:backward/:auto
    end
    struct Absorbing <: HomogeneousBoundaryCondition end