# General dispatcher
"""
   `diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition)`

Return the diffusion operators `(L_1_minus, L_1_plus, L_2)` w.r.t the supplied grid and BCs.

`x` is a grid (either an `AbstractRange`, in which we use specialized uniform grid code, or an `AbstractArray`). The
first BC binds at the lower end of the grid (i.e., `x[1]`), and the latter at the high end. The BCs are either a `Reflecting()`,
or "Dirichlet" boundary condition `v'(x) = 0`, or `Mixed(x::T) where {T <: Real}`, corresponding to "Robin" boundary conditions.
"""
diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition) = _diffusionoperators(x, BC1, BC2)

# Uniform grid, (Reflecting, Reflecting)
function _diffusionoperators(x::AbstractRange, BC1::Reflecting, BC2::Reflecting)
    T = eltype(x) # get data type of the range
    Δ = step(x)
    P = length(x)

    # define L_1_minus
    dl_m1 = -ones(T, P-1)
    d_m1 = ones(T, P)
    du_m1 = zeros(T, P-1)
    d_m1[1] = d_m1[1] + dl_m1[1]
    L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)/Δ

    # define L_1_plus
    dl_1 = zeros(T, P-1)
    d_1 = -ones(T, P)
    du_1 = ones(T, P-1)
    d_1[end] = d_1[end] + du_1[end]
    L_1_plus = Tridiagonal(dl_1, d_1, du_1)/Δ

    # define L_2
    dl_2 = ones(T, P-1)
    d_2 = -2 * ones(T, P)
    d_2[1] = -one(T)
    d_2[end] = -one(T)
    du_2 = ones(T, P-1)
    L_2 = Tridiagonal(dl_2, d_2, du_2)/(Δ^2)

    # return
    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
end

# Uniform grid, (Mixed, Mixed)
function _diffusionoperators(x::AbstractRange, BC1::Mixed, BC2::Mixed)
    Δ = step(x)
    ξ_lb = BC1.ξ
    ξ_ub = BC2.ξ

    # extract diffusion operators with reflecting barrier conditions first
    L_1_minus, L_1_plus, L_2 = diffusionoperators(x, Reflecting(), Reflecting())

    # apply boundary condition constraints
    L_1_minus[1,1] -= ξ_lb * Δ
    L_1_plus[end,end] -= ξ_ub * Δ
    L_2[1,1] += ξ_lb * Δ
    L_2[end,end] -= ξ_ub * Δ

    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
end

# Irregular grid, (Reflecting, Reflecting)
function _diffusionoperators(x::AbstractArray, BC1::Reflecting, BC2::Reflecting)
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
   d_m1[1] = d_m1[1] + dl_m1[1]
   du_m1 = zeros(T, P-1)./Δ_p[1:end-1]
   L_1_minus = Tridiagonal(dl_m1, d_m1, du_m1)

   # define L_1_plus
   dl_p1 = zeros(T, P-1)./Δ_m[2:end]
   d_p1 = -ones(T, P)./Δ_p
   du_p1 = ones(T, P-1)./Δ_p[1:end-1]
   d_p1[end] = d_p1[end] + du_p1[end]
   L_1_plus = Tridiagonal(dl_p1, d_p1, du_p1)

   # define L_2
   Δ=Δ_p+Δ_m
   dl_2 = 2*ones(T, P-1)./(Δ_m[2:end].*Δ[2:end])
   d_2 = -2*ones(T, P)
   d_2[1] = -one(T)
   d_2[end] = -one(T)
   d_2 = d_2./(Δ_p.*Δ_m)
   du_2 = 2*ones(T, P-1)./(Δ_p[1:end-1].*Δ[1:end-1])
   L_2 = Tridiagonal(dl_2, d_2, du_2)

   # return
   return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
end

# Irregular grid, (Mixed, Mixed)
function _diffusionoperators(x::AbstractArray, BC1::Mixed, BC2::Mixed)
    d = diff(x)
    Δ_1 = d[1]
    Δ_P = d[end]
    ξ_lb = BC1.ξ
    ξ_ub = BC2.ξ

    # extract diffusion operators with reflecting barrier conditions first
    L_1_minus, L_1_plus, L_2 = diffusionoperators(x, Reflecting(), Reflecting())

    # apply boundary condition constraints
    L_1_minus[1,1] -= ξ_lb * Δ_1
    L_1_plus[end,end] -= ξ_ub * Δ_P
    L_2[1,1] += ξ_lb * Δ_1
    L_2[end,end] -= ξ_ub * Δ_P

    return (L_1_minus = L_1_minus, L_1_plus = L_1_plus, L_2 = L_2)
end
