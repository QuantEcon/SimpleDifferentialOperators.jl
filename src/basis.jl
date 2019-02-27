# returns matrices that define the basis for interiors of operators for 
# L_1 based on BD, L_1 based on FD, and L_2 based on CD respectively
# where P is the number of nodes
function get_basis_L₁₋(x::AbstractRange)
    T = eltype(x) # get data type of the range
    Δ = step(x)
    M = length(x)

    # define L₁₋
    dl_m1 = -ones(T, M-1)
    d_m1 = ones(T, M)
    du_m1 = zeros(T, M-1)
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)/Δ
end

function get_basis_L₁₊(x::AbstractRange)
    T = eltype(x) # get data type of the range
    Δ = step(x)
    M = length(x)

    # define L₁₊
    dl_1 = zeros(T, M-1)
    d_1 = -ones(T, M)
    du_1 = ones(T, M-1)
    L₁₊ = Tridiagonal(dl_1, d_1, du_1)/Δ
end

function get_basis_L₂(x::AbstractRange)
    T = eltype(x) # get data type of the range
    Δ = step(x)
    M = length(x)

    # define L₂
    dl_2 = ones(T, M-1)
    d_2 = -2 * ones(T, M)
    du_2 = ones(T, M-1)
    L₂ = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)
end

function get_basis_L₁₋(x::AbstractArray)
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    # define L₁₋
    dl_m1 = -ones(T, M-1)./Δ_m[2:end]
    d_m1 = ones(T, M)./Δ_m
    du_m1 = zeros(T, M-1)./Δ_p[1:end-1]
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)
end

function get_basis_L₁₊(x::AbstractArray)
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    # define L₁₊
    dl_p1 = zeros(T, M-1)./Δ_m[2:end]
    d_p1 = -ones(T, M)./Δ_p
    du_p1 = ones(T, M-1)./Δ_p[1:end-1]
    L₁₊ = Tridiagonal(dl_p1, d_p1, du_p1)
 
end

function get_basis_L₂(x::AbstractArray)
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d
    # define L₂
    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(T, M-1)./(Δ_m[2:end].*Δ[2:end])
    d_2 = -2*ones(T, M)
    d_2 = d_2./(Δ_p.*Δ_m)
    du_2 = 2*ones(T, M-1)./(Δ_p[1:end-1].*Δ[1:end-1])
    L₂ = Tridiagonal(dl_2, d_2, du_2)
end

# convenient functions that returns basis corresponding to discretization method given by `method`
get_basis_operator(x, method::BackwardFirstDifference) = get_basis_L₁₋(x)
get_basis_operator(x, method::ForwardFirstDifference) = get_basis_L₁₊(x)
get_basis_operator(x, method::CentralSecondDifference) = get_basis_L₂(x)