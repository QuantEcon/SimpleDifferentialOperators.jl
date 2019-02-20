# returns a tuple of P by P matrices `(L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)`
# that defines the basis for interiors of operators for 
# L_1 based on BD, L_1 based on FD, and L_2 based on CD respectively
# where P is the number of nodes
function get_operator_basis(x::AbstractRange)
    T = eltype(x) # get data type of the range
    Δ = step(x)
    P = length(x)

    # define L_1_minus
    dl_m1 = -ones(T, P-1)
    d_m1 = ones(T, P)
    du_m1 = zeros(T, P-1)
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)/Δ

    # define L_1_plus
    dl_1 = zeros(T, P-1)
    d_1 = -ones(T, P)
    du_1 = ones(T, P-1)
    L_1_plus = Tridiagonal(dl_1, d_1, du_1)/Δ

    # define L_2
    dl_2 = ones(T, P-1)
    d_2 = -2 * ones(T, P)
    du_2 = ones(T, P-1)
    L_2 = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)

    # return
    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
end

# basis for irregular grids
function get_operator_basis(x::AbstractArray)
    # define preliminaries
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    P = length(x)
    Δ_m = zeros(T, P)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, P)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d
 
    # define L_1_minus
    dl_m1 = -ones(T, P-1)./Δ_m[2:end]
    d_m1 = ones(T, P)./Δ_m
    du_m1 = zeros(T, P-1)./Δ_p[1:end-1]
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)
 
    # define L_1_plus
    dl_p1 = zeros(T, P-1)./Δ_m[2:end]
    d_p1 = -ones(T, P)./Δ_p
    du_p1 = ones(T, P-1)./Δ_p[1:end-1]
    L_1_plus = Tridiagonal(dl_p1, d_p1, du_p1)
 
    # define L_2
    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(T, P-1)./(Δ_m[2:end].*Δ[2:end])
    d_2 = -2*ones(T, P)
    d_2 = d_2./(Δ_p.*Δ_m)
    du_2 = 2*ones(T, P-1)./(Δ_p[1:end-1].*Δ[1:end-1])
    L_2 = Tridiagonal(dl_2, d_2, du_2)
 
    # return
    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
 end