# Concrete "under the hood" methods.
# (NoBoundary, NoBoundary)
function DifferentialOperator(x::AbstractRange, bc::Tuple{NoBoundary, NoBoundary}, method::ForwardFirstDifference)
    T = eltype(x)
    d = diff(x)
    Δ = step(x)
    M = length(x)

    dl_1 = zeros(T, M-1)
    d_1 = -ones(T, M)
    du_1 = ones(T, M-1)
    L₁₊ = Tridiagonal(dl_1, d_1, du_1)/Δ

    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)
    L₁₊ = sparse([col_zeros L₁₊ col_ub/Δ])
end

function DifferentialOperator(x::AbstractRange, bc::Tuple{NoBoundary, NoBoundary}, method::BackwardFirstDifference)
    T = eltype(x)
    d = diff(x)
    Δ = step(x)
    M = length(x)

    dl_m1 = -ones(T, M-1)
    d_m1 = ones(T, M)
    du_m1 = zeros(T, M-1)
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)/Δ

    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)
    L₁₋ = sparse([-col_lb/Δ L₁₋ col_zeros])
end

function DifferentialOperator(x::AbstractRange, bc::Tuple{NoBoundary, NoBoundary}, method::CentralSecondDifference)
    T = eltype(x)
    d = diff(x)
    Δ = step(x)
    M = length(x)

    dl_2 = ones(T, M-1)
    d_2 = -2 * ones(T, M)
    du_2 = ones(T, M-1)
    L₂ = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)

    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)
    L₂ = sparse([col_lb/(Δ*Δ) L₂ col_ub/(Δ*Δ)])
end

function DifferentialOperator(x, bc::Tuple{NoBoundary, NoBoundary}, method::ForwardFirstDifference)
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    M = length(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    # define L_1_minus
    dl_p1 = zeros(T, M-1)./Δ_m[2:end]
    d_p1 = -ones(T, M)./Δ_p
    du_p1 = ones(T, M-1)./Δ_p[1:end-1]
    L₁₊ = Tridiagonal(dl_p1, d_p1, du_p1)

    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)
    L₁₊ = sparse([col_zeros L₁₊ col_ub/Δ_M])
end

function DifferentialOperator(x, bc::Tuple{NoBoundary, NoBoundary}, method::BackwardFirstDifference)
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    M = length(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    # define L_1_minus
    dl_m1 = -ones(T, M-1)./Δ_m[2:end]
    d_m1 = ones(T, M)./Δ_m
    du_m1 = zeros(T, M-1)./Δ_p[1:end-1]
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)

    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)
    L₁₋ = sparse([-col_lb/Δ_1 L₁₋ col_zeros])
end

function DifferentialOperator(x, bc::Tuple{NoBoundary, NoBoundary}, method::CentralSecondDifference)
    T = eltype(x) # get data type of the grid
    d = diff(x) # using the first difference as diff from ghost node
    Δ_1 = d[1]
    Δ_M = d[end]
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(T, M-1)./(Δ_m[2:end].*Δ[2:end])
    d_2 = -2*ones(T, M)
    d_2 = d_2./(Δ_p.*Δ_m)
    du_2 = 2*ones(T, M-1)./(Δ_p[1:end-1].*Δ[1:end-1])
    L₂ = Tridiagonal(dl_2, d_2, du_2)

    col_zeros = zeros(T, M)
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = one(T)
    col_ub[end] = one(T)
    L₂ = sparse([col_lb/(Δ_1*Δ_1) L₂ col_ub/(Δ_M*Δ_M)])
end


# (Reflecting, Reflecting)
function DifferentialOperator(x::AbstractRange, bc::Tuple{Reflecting, Reflecting}, method::ForwardFirstDifference)
    T = eltype(x)
    Δ = step(x)
    M = length(x)

    dl_1 = zeros(T, M-1)
    d_1 = -ones(T, M)
    du_1 = ones(T, M-1)
    L₁₊ = Tridiagonal(dl_1, d_1, du_1)/Δ
    L₁₊[end, end] = zero(T)
    return L₁₊
end

function DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::ForwardFirstDifference)
    T = eltype(x)
    d = diff(x)
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    dl_p1 = zeros(T, M-1)./Δ_m[2:end]
    d_p1 = -ones(T, M)./Δ_p
    du_p1 = ones(T, M-1)./Δ_p[1:end-1]
    L₁₊ = Tridiagonal(dl_p1, d_p1, du_p1)
    L₁₊[end, end] = zero(T)
    return L₁₊
end

function DifferentialOperator(x::AbstractRange, bc::Tuple{Reflecting, Reflecting}, method::BackwardFirstDifference)
    T = eltype(x)
    Δ = step(x)
    M = length(x)

    dl_m1 = -ones(T, M-1)
    d_m1 = ones(T, M)
    du_m1 = zeros(T, M-1)
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)/Δ
    L₁₋[1, 1] = zero(T)
    return L₁₋
end

function DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::BackwardFirstDifference)
    T = eltype(x)
    d = diff(x)
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    dl_m1 = -ones(T, M-1)./Δ_m[2:end]
    d_m1 = ones(T, M)./Δ_m
    du_m1 = zeros(T, M-1)./Δ_p[1:end-1]
    L₁₋ = Tridiagonal(dl_m1, d_m1, du_m1)
    L₁₋[1, 1] = zero(T)
    return L₁₋
end

function DifferentialOperator(x::AbstractRange, bc::Tuple{Reflecting, Reflecting}, method::CentralSecondDifference)
    T = eltype(x)
    Δ = step(x)
    M = length(x)

    dl_2 = ones(T, M-1)
    d_2 = -2 * ones(T, M)
    du_2 = ones(T, M-1)
    L₂ = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)
    L₂[1,1] /= 2
    L₂[end,end] /= 2
    return L₂
end

function DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::CentralSecondDifference)
    T = eltype(x)
    d = diff(x)
    M = length(x)
    Δ_m = zeros(T, M)
    Δ_m[1] = d[1]
    Δ_m[2:end] = d
    Δ_p = zeros(T, M)
    Δ_p[end] = d[end]
    Δ_p[1:end-1] = d

    Δ=Δ_p+Δ_m
    dl_2 = 2*ones(T, M-1)./(Δ_m[2:end].*Δ[2:end])
    d_2 = -2*ones(T, M)
    d_2 = d_2./(Δ_p.*Δ_m)
    du_2 = 2*ones(T, M-1)./(Δ_p[1:end-1].*Δ[1:end-1])
    L₂ = Tridiagonal(dl_2, d_2, du_2)
    L₂[1,1] /= 2
    L₂[end,end] /= 2
    return L₂
end

# (Mixed, Mixed)
function DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::ForwardFirstDifference)
    d = diff(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ
    L₁₊ = DifferentialOperator(x, (Reflecting(), Reflecting()), method)
    L₁₊[end, end] -= ξ_ub
    return L₁₊
end

function DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::BackwardFirstDifference)
    d = diff(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ
    L₁₋ = DifferentialOperator(x, (Reflecting(), Reflecting()), method)
    L₁₋[1, 1] -= ξ_lb
    return L₁₋
end

function DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::CentralSecondDifference)
    d = diff(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ
    L₂ = DifferentialOperator(x, (Reflecting(), Reflecting()), method)
    L₂[1, 1] += ξ_lb / Δ_1
    L₂[end, end] -= ξ_ub / Δ_M
    return L₂
end

# Convenience calls
    L₁₋(x, bc) = DifferentialOperator(x, bc, BackwardFirstDifference())
    L₁₊(x, bc) = DifferentialOperator(x, bc, ForwardFirstDifference())
    L₂(x, bc) = DifferentialOperator(x, bc, CentralSecondDifference())

    L̄₁₋(x) = DifferentialOperator(x, (NoBoundary(), NoBoundary()), BackwardFirstDifference())
    L̄₁₊(x) = DifferentialOperator(x, (NoBoundary(), NoBoundary()), ForwardFirstDifference())
    L̄₂(x)  = DifferentialOperator(x, (NoBoundary(), NoBoundary()), CentralSecondDifference())

    function x̄(x)
        d = diff(x) # dispatches based on AbstractArray or not
        x̄ = collect([x[1] - d[1]; x; x[end] + d[end]])
    end

    diffusionoperators(x, bc) = (L₁₋ = L₁₋(x, bc), L₁₊ = L₁₊(x, bc), L₂ = L₂(x, bc), x̄ = x̄(x))
