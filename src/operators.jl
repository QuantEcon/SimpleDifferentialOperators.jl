# General dispatcher
"""
   `diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition)`

Returns a tuple of diffusion operators and extended grid `(L_1_minus, L_1_plus, L_2, x_bar)`
with specified boundary conditions.

Given a grid `x` of length `M`, return diffusion operators for negative drift, positive drift,
and central differences. BC1 is applied to the lower bound, and BC2 to the upper. `x_bar` is a `(M+2)` array that
represents the extended grid whose first and last elements represent the ghost nodes
just before `x[1]` and `x[end]`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())
(L_1_minus = [0.0 0.0 0.0; -1.0 1.0 0.0; 0.0 -1.0 1.0], L_1_plus = [-1.0 1.0 0.0; 0.0 -1.0 1.0; 0.0 0.0 0.0], L_2 = [-1.0 1.0 0.0; 1.0 -2.0 1.0; 0.0 1.0 -1.0], x_bar = [0, 1, 2, 3, 4])
```
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
    Δ_M = d[end]
    ξ_lb = BC1.ξ
    ξ_ub = BC2.ξ

    # extract diffusion operators with reflecting barrier conditions first
    L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())

    # apply boundary condition constraints
    L_1_minus[1,1] -= ξ_lb
    L_1_plus[end,end] -= ξ_ub
    L_2[1,1] += ξ_lb / Δ_1
    L_2[end,end] -= ξ_ub / Δ_M

    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
end
