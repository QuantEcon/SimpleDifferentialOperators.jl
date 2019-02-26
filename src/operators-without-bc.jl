"""
   `diffusionoperators(x)`

Returns a tuple of diffusion operators and extended grid `(L_1_minus, L_1_plus, L_2, x_bar)`
without applying any boundary conditions.

Given a grid `x` of length `M`, return diffusion operators for negative drift, positive drift,
and central differences. No boundary conditions are applied. `x_bar` is a `(M+2)` array that
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
diffusionoperators(x) = diffusionoperators(x, NoBoundary())

function diffusionoperators(x, BC::NoBoundary)
    T = eltype(x) # get data type of the range
    d = diff(x)
    M = length(x)

    # get grid diff to be used for ghost nodes
    Δ_1 = d[1]
    Δ_M = d[end]
    L_1_minus, L_1_plus, L_2, x_bar = get_operator_basis(x)

    # get columns for ghost nodes on lower bound and upper bound
    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)

    # attach the corresponding columns to operator basis matrices
    L_1_minus = sparse([-col_lb/Δ_1 L_1_minus col_zeros])
    L_1_plus = sparse([col_zeros L_1_plus col_ub/Δ_M])
    L_2 = sparse([col_lb/(Δ_1*Δ_1) L_2 col_ub/(Δ_M*Δ_M)])

    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
end
