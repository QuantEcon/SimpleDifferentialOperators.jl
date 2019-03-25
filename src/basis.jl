# returns matrices that define the basis for interiors of operators for 
# L_1 based on BD, L_1 based on FD, and L_2 based on CD respectively
# where P is the number of nodes
function get_basis_L₁₋(x̄::AbstractRange)
    T = eltype(x̄) # get data type of the range
    Δ = step(x̄)
    M = length(x̄) - 2

    # define L₁₋
    dl_m1 = -ones(T, M-1)
    d_m1 = ones(T, M)
    du_m1 = zeros(T, M-1)
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)/Δ
end

function get_basis_L₁₊(x̄::AbstractRange)
    T = eltype(x̄) # get data type of the range
    Δ = step(x̄)
    M = length(x̄) - 2

    # define L₁₊
    dl_1 = zeros(T, M-1)
    d_1 = -ones(T, M)
    du_1 = ones(T, M-1)
    L₁₊ = Tridiagonal(dl_1, d_1, du_1)/Δ
end

function get_basis_L₂(x̄::AbstractRange)
    T = eltype(x̄) # get data type of the range
    Δ = step(x̄)
    M = length(x̄) - 2

    # define L₂
    dl_2 = ones(T, M-1)
    d_2 = -2 * ones(T, M)
    du_2 = ones(T, M-1)
    L₂ = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)
end

function get_basis_L₁₋(x̄::AbstractArray)
    T = eltype(x̄) # get data type of the grid
    d = diff(x̄) # using the first difference as diff from ghost node
    M = length(x̄) - 2
    Δ₋ = d[1:end-1]
    Δ₊ = d[2:end]

    # define L₁₋
    dl_m1 = -ones(T, M-1)./Δ₋[2:end]
    d_m1 = ones(T, M)./Δ₋
    du_m1 = zeros(T, M-1)./Δ₊[1:end-1]
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)
end

function get_basis_L₁₊(x̄::AbstractArray)
    T = eltype(x̄) # get data type of the grid
    d = diff(x̄) # using the first difference as diff from ghost node
    M = length(x̄) - 2
    Δ₋ = d[1:end-1]
    Δ₊ = d[2:end]

    # define L₁₊
    dl_p1 = zeros(T, M-1)./Δ₋[2:end]
    d_p1 = -ones(T, M)./Δ₊
    du_p1 = ones(T, M-1)./Δ₊[1:end-1]
    L₁₊ = Tridiagonal(dl_p1, d_p1, du_p1)
end

function get_basis_L₂(x̄::AbstractArray)
    T = eltype(x̄) # get data type of the grid
    d = diff(x̄) # using the first difference as diff from ghost node
    M = length(x̄) - 2
    Δ₋ = d[1:end-1]
    Δ₊ = d[2:end]
    
    # define L₂
    Δ=Δ₊+Δ₋
    dl_2 = 2*ones(T, M-1)./(Δ₋[2:end].*Δ[2:end])
    d_2 = -2*ones(T, M)./(Δ₊.*Δ₋)
    du_2 = 2*ones(T, M-1)./(Δ₊[1:end-1].*Δ[1:end-1])
    L₂ = Tridiagonal(dl_2, d_2, du_2)
end

# convenience functions that returns basis corresponding to discretization method given by `method`
get_basis_operator(x̄, method::BackwardFirstDifference) = get_basis_L₁₋(x̄)
get_basis_operator(x̄, method::ForwardFirstDifference) = get_basis_L₁₊(x̄)
get_basis_operator(x̄, method::CentralSecondDifference) = get_basis_L₂(x̄)