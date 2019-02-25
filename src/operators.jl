# General dispatcher
"""
   `diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition)`

Return the diffusion operators `(L_1_minus, L_1_plus, L_2)` w.r.t the supplied grid and BCs.

Given a grid `x` of length `M`, return `L_1_minus`, `L_1_plus`, `L_2` that are M by M matrices representing 
L_1 based on BD, L_1 based on FD, and L_2 based on CD respectively, 
where a lower boundary condition `BC1` and upper boundary condition `BC2` are applied.
`x_bar` is a `(M+2)` array that represents the extended grid whose first element
and the last element represent the ghost nodes on lower boundary and upper boundary.

# Examples
```jldoctest
julia> x = 1:3
1:3

julia> L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x)
(L_1_minus =
  [1, 1]  =  -1.0
  [1, 2]  =  1.0
  [2, 2]  =  -1.0
  [1, 3]  =  0.0
  [2, 3]  =  1.0
  [3, 3]  =  -1.0
  [2, 4]  =  0.0
  [3, 4]  =  1.0, L_1_plus =
  [1, 2]  =  -1.0
  [2, 2]  =  0.0
  [1, 3]  =  1.0
  [2, 3]  =  -1.0
  [3, 3]  =  0.0
  [2, 4]  =  1.0
  [3, 4]  =  -1.0
  [3, 5]  =  1.0, L_2 =
  [1, 1]  =  1.0
  [1, 2]  =  -2.0
  [2, 2]  =  1.0
  [1, 3]  =  1.0
  [2, 3]  =  -2.0
  [3, 3]  =  1.0
  [2, 4]  =  1.0
  [3, 4]  =  -2.0
  [3, 5]  =  1.0, x_bar = [0, 1, 2, 3, 4])
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
