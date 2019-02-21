"""
   `diffusionoperators(x)`

Returns a tuple of diffusion operators and extended grid `(L_1_minus, L_1_plus, L_2, x_bar)` 
where `L_1_minus`, `L_1_plus`, `L_2` are P by (P+2) matrices that represent 
L_1 based on BD, L_1 based on FD, and L_2 based on CD respectively, without any boundary conditions
where P is `length(x)`. `x_bar` is `(P+2)` array that represents the extended grid whose first element
and the last element represent the ghost nodes on lower boundary and upper boundary, respectively.
"""
diffusionoperators(x) = diffusionoperators(x, NoBoundary())

function diffusionoperators(x, BC::NoBoundary)
    T = eltype(x) # get data type of the range
    d = diff(x)
    P = length(x)

    # get grid diff to be used for ghost nodes
    Δ_1 = d[1]
    Δ_P = d[end]
    L_1_minus, L_1_plus, L_2 = get_operator_basis(x)

    # get columns for ghost nodes on lower bound and upper bound
    col_zeros = zeros(T, P)
    col_lb = zeros(T, P)
    col_ub = zeros(T, P)
    col_lb[1] = one(T)
    col_ub[end] = one(T)

    # attach the corresponding columns to operator basis matrices
    L_1_minus = sparse([-col_lb/Δ_1 L_1_minus col_zeros])
    L_1_plus = sparse([col_zeros L_1_plus col_ub/Δ_P])
    L_2 = sparse([col_lb/(Δ_1*Δ_1) L_2 col_ub]/(Δ_P*Δ_P))

    # define extended grid with ghost nodes
    x_bar = collect([x[1] - Δ_1; x; x[end] + Δ_P])
    
    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2, x_bar = x_bar)
end