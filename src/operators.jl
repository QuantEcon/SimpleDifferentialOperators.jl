# General dispatcher
"""
   `diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition)`

Return the diffusion operators `(L_1_minus, L_1_plus, L_2)` w.r.t the supplied grid and BCs.

`x` is a grid (either an `AbstractRange`, in which we use specialized uniform grid code, or an `AbstractArray`). The
first BC binds at the lower end of the grid (i.e., `x[1]`), and the latter at the high end. The BCs are either a `Reflecting()`,
or "Dirichlet" boundary condition `v'(x) = 0`, or `Mixed(x::T) where {T <: Real}`, corresponding to "Robin" boundary conditions.
"""
diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition) = _diffusionoperators(x, BC1, BC2)

# (Reflecting, Reflecting)
function _diffusionoperators(x, BC1::Reflecting, BC2::Reflecting)
    T = eltype(x)
    L_1_minus, L_1_plus, L_2, x_bar = get_operator_basis(x)

    # apply boundary conditions
    L_1_minus[1,1] = zero(T)
    L_1_plus[end,end] = zero(T)
    L_2[1,1] /= 2
    L_2[end,end] /= 2 

    # return
    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
end


# (Mixed, Mixed)
function _diffusionoperators(x, BC1::Mixed, BC2::Mixed)
    d = diff(x)
    Δ_1 = d[1]
    Δ_P = d[end]
    ξ_lb = BC1.ξ
    ξ_ub = BC2.ξ

    # extract diffusion operators with reflecting barrier conditions first
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())

    # apply boundary condition constraints
    L_1_minus[1,1] -= ξ_lb
    L_1_plus[end,end] -= ξ_ub
    L_2[1,1] += ξ_lb / Δ_1
    L_2[end,end] -= ξ_ub / Δ_P

    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
end
