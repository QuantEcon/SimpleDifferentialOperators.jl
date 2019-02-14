# General dispatcher
diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition) = _diffusionoperators(x, BC1, BC2)

#=
    UNIFORM GRID
=#

# Uniform grid, two reflecting
function _diffusionoperators(x::AbstractRange, BC1::Reflecting, BC2::Reflecting)
    return (x = "Uniform", BC1 = "Reflecting", BC2 = "Reflecting")
end

# Uniform grid, (Reflecting, Mixed)
function _diffusionoperators(x::AbstractRange, BC1::Reflecting, BC2::Mixed)
    return (x = "Uniform", BC1 = "Reflecting", BC2 = "Mixed")
end

# Uniform grid, (Mixed, Reflecting)
function _diffusionoperators(x::AbstractRange, BC1::Mixed, BC2::Reflecting)
    return (x = "Uniform", BC1 = "Mixed", BC2 = "Reflecting")
end

# Uniform grid, (Mixed, Mixed)
function _diffusionoperators(x::AbstractRange, BC1::Mixed, BC2::Mixed)
    return (x = "Uniform", BC1 = "Mixed", BC2 = "Mixed")
end

# Uniform grid, (Reflecting, Absorbing)
function _diffusionoperators(x::AbstractRange, BC1::Reflecting, BC2::Absorbing)
    return (x = "Uniform", BC1 = "Reflecting", BC2 = "Absorbing")
end

# Uniform grid, (Absorbing, Reflecting)
function _diffusionoperators(x::AbstractRange, BC1::Absorbing, BC2::Reflecting)
    return (x = "Uniform", BC1 = "Absorbing", BC2 = "Reflecting")
end

# Uniform grid, (Absorbing, Absorbing)
function _diffusionoperators(x::AbstractRange, BC1::Absorbing, BC2::Absorbing)
    return (x = "Uniform", BC1 = "Absorbing", BC2 = "Absorbing")
end

# Uniform grid, (Absorbing, Mixed)
function _diffusionoperators(x::AbstractRange, BC1::Absorbing, BC2::Mixed)
    return (x = "Uniform", BC1 = "Absorbing", BC2 = "Mixed")
end

# Uniform grid, (Mixed, Absorbing)
function _diffusionoperators(x::AbstractRange, BC1::Mixed, BC2::Absorbing)
    return (x = "Uniform", BC1 = "Mixed", BC2 = "Absorbing")
end

#=
    IRREGULAR GRID
=#

# Irregular grid, two reflecting
function _diffusionoperators(x::AbstractArray, BC1::Reflecting, BC2::Reflecting)
    return (x = "Irregular", BC1 = "Reflecting", BC2 = "Reflecting")
end

# Irregular grid, (Reflecting, Mixed)
function _diffusionoperators(x::AbstractArray, BC1::Reflecting, BC2::Mixed)
    return (x = "Irregular", BC1 = "Reflecting", BC2 = "Mixed")
end

# Irregular grid, (Mixed, Reflecting)
function _diffusionoperators(x::AbstractArray, BC1::Mixed, BC2::Reflecting)
    return (x = "Irregular", BC1 = "Mixed", BC2 = "Reflecting")
end

# Irregular grid, (Mixed, Mixed)
function _diffusionoperators(x::AbstractArray, BC1::Mixed, BC2::Mixed)
    return (x = "Irregular", BC1 = "Mixed", BC2 = "Mixed")
end

# Irregular grid, (Reflecting, Absorbing)
function _diffusionoperators(x::AbstractArray, BC1::Reflecting, BC2::Absorbing)
    return (x = "Irregular", BC1 = "Reflecting", BC2 = "Absorbing")
end

# Irregular grid, (Absorbing, Reflecting)
function _diffusionoperators(x::AbstractArray, BC1::Absorbing, BC2::Reflecting)
    return (x = "Irregular", BC1 = "Absorbing", BC2 = "Reflecting")
end

# Irregular grid, (Absorbing, Absorbing)
function _diffusionoperators(x::AbstractArray, BC1::Absorbing, BC2::Absorbing)
    return (x = "Irregular", BC1 = "Absorbing", BC2 = "Absorbing")
end

# Irregular grid, (Absorbing, Mixed)
function _diffusionoperators(x::AbstractArray, BC1::Absorbing, BC2::Mixed)
    return (x = "Irregular", BC1 = "Absorbing", BC2 = "Mixed")
end

# Irregular grid, (Mixed, Absorbing)
function _diffusionoperators(x::AbstractArray, BC1::Mixed, BC2::Absorbing)
    return (x = "Irregular", BC1 = "Mixed", BC2 = "Absorbing")
end
