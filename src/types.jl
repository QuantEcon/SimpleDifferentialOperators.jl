# Abstract types.
abstract type BoundaryCondition end
abstract type HomogeneousBoundaryCondition <: BoundaryCondition end
abstract type InhomogeneousBoundaryCondition <: BoundaryCondition end

# Concrete types.
struct Reflecting <: HomogeneousBoundaryCondition end
struct Mixed{T} <: HomogeneousBoundaryCondition where {T <: Real}
    ξ::T
end
struct Absorbing{T1, T2} <: InhomogeneousBoundaryCondition where {T1 <: Real, T2 <: Real}
    x::T1
    y::T2
end
struct NoBoundary <: BoundaryCondition end