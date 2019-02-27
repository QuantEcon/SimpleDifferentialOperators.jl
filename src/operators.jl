# Concrete "under the hood" methods.
# (NoBoundary, NoBoundary)
function DifferentialOperator(x, bc::Tuple{NoBoundary, NoBoundary}, method::DifferenceMethod)
    T = eltype(x)
    d = diff(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    M = length(x)

    # get basis operator on interior nodes
    L_basis = get_basis_operator(x, method)

    # add columns for ghost nodes next to boundaries
    col_lb = zeros(T, M)
    col_ub = zeros(T, M)
    col_lb[1] = typeof(method) <: BackwardFirstDifference ? -(one(T) / Δ_1) : zero(T)
    col_lb[1] = typeof(method) <: CentralSecondDifference ? (one(T) / (Δ_1*Δ_1)) : col_lb[1]
    col_ub[end] = typeof(method) <: ForwardFirstDifference ? (one(T) / Δ_M) : zero(T)
    col_ub[end] = typeof(method) <: CentralSecondDifference ? (one(T) / (Δ_M*Δ_M)) : col_ub[end]

    L = sparse([col_lb L_basis col_ub])
end

# (Reflecting, Reflecting)
function DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)
    T = eltype(x)

    # get basis operator on interior nodes
    L = get_basis_operator(x, method)

    # apply boundary conditions
    L[1,1] = typeof(method) <: BackwardFirstDifference ? zero(T) : L[1,1]
    L[1,1] = typeof(method) <: CentralSecondDifference ? (L[1,1] / 2) : L[1,1]
    L[end,end] = typeof(method) <: ForwardFirstDifference ? zero(T) : L[end,end]
    L[end,end] = typeof(method) <: CentralSecondDifference ? (L[end,end] / 2) : L[end,end] 

    return L
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
