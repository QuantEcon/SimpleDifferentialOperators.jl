# General dispatcher
diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition) = _diffusionoperators(x, BC1, BC2)

# Uniform grid, (Reflecting, Reflecting)
function _diffusionoperators(x::AbstractRange, BC1::Reflecting, BC2::Reflecting)
    return (x = "Uniform", BC1 = "Reflecting", BC2 = "Reflecting")
end

# Uniform grid, (Mixed, Mixed)
function _diffusionoperators(x::AbstractRange, BC1::Mixed, BC2::Mixed)
    return (x = "Uniform", BC1 = "Mixed", BC2 = "Mixed")
end

# Irregular grid, (Reflecting, Reflecting)
function _diffusionoperators(x::AbstractArray, BC1::Reflecting, BC2::Reflecting)
    return (x = "Irregular", BC1 = "Reflecting", BC2 = "Reflecting")
end

# Irregular grid, (Mixed, Mixed)
function _diffusionoperators(x::AbstractArray, BC1::Mixed, BC2::Mixed)
    return (x = "Irregular", BC1 = "Mixed", BC2 = "Mixed")
end
