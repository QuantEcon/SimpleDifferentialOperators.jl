function L₀affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    return zeros(M)
end

function L₁₋affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    b = zeros(M)

    # apply affine operator if needed
    Δ_1m = x̄[2] - x̄[1]
    b[1] = (typeof(bc[1]) <: NonhomogeneousAbsorbing) ? bc[1].S / Δ_1m : b[1]
    return b
end

function L₁₊affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    b = zeros(M)

    # apply affine operator if needed
    Δ_Mp = x̄[end] - x̄[end-1]
    b[end] = (typeof(bc[2]) <: NonhomogeneousAbsorbing) ? -bc[2].S / Δ_Mp : b[end]
    return b
end

function L₂affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    b = zeros(M)

    # apply affine operator if needed
    Δ_1p = x̄[3] - x̄[2]
    Δ_1m = x̄[2] - x̄[1]
    Δ_Mp = x̄[end] - x̄[end-1]
    Δ_Mm = x̄[end-1] - x̄[end-2]
    b[1] = (typeof(bc[1]) <: NonhomogeneousAbsorbing) ? -2*bc[1].S/(Δ_1m*(Δ_1m+Δ_1p)) : b[1]
    b[end] = (typeof(bc[2]) <: NonhomogeneousAbsorbing) ? -2*bc[2].S/(Δ_Mp*(Δ_Mm+Δ_Mp)) : b[end]

    return b
end

# take M by (M+2) extended operator 
function Laffine(L, bc::Tuple{BoundaryCondition, BoundaryCondition})
    M = size(L)[1]
    b = zeros(M)
    
    # apply affine operator if needed
    b[1] = (typeof(bc[1]) <: NonhomogeneousAbsorbing) ? -sum(L[:,1])*bc[1] : b[1]
    b[end] = (typeof(bc[2]) <: NonhomogeneousAbsorbing) ? -sum(L[:,end])*bc[end] : b[end]

    return b
end