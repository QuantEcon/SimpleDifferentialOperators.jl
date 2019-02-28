"""
    `ExtensionDifferentialOperator(x, method::DifferenceMethod)`
Returns a discretized differential operator of `length(x)` by `length(x) + 2` matrix
whose first and last columns are applied to the ghost nodes just before `x[1]` and `x[end]` respectively
under no boundary condition using finite difference method specified by `method`.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> ExtensionDifferentialOperator(x, BackwardFirstDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
  0.0   0.0   ⋅
 -1.0   1.0  0.0
   ⋅   -1.0  1.0

julia> ExtensionDifferentialOperator(x, ForwardFirstDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0   ⋅
  0.0  -1.0  1.0
   ⋅    0.0  0.0

julia> ExtensionDifferentialOperator(x, CentralSecondDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅
  1.0  -2.0   1.0
   ⋅    1.0  -1.0
```
"""
function ExtensionDifferentialOperator(x, method::DifferenceMethod)
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

"""
    `DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)`
Returns a discretized differential operator of `length(x)` by `length(x)` matrix
under reflecting boundary conditions from `bc` using finite difference method specified by `method`. 
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> DifferentialOperator(x, (Reflecting(), Reflecting()), BackwardFirstDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
  0.0   0.0   ⋅
 -1.0   1.0  0.0
   ⋅   -1.0  1.0

julia> DifferentialOperator(x, (Reflecting(), Reflecting()), ForwardFirstDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0   ⋅
  0.0  -1.0  1.0
   ⋅    0.0  0.0

julia> DifferentialOperator(x, (Reflecting(), Reflecting()), CentralSecondDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅
  1.0  -2.0   1.0
   ⋅    1.0  -1.0
```
"""
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

"""
    `DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)`
Returns a discretized differential operator of `length(x)` by `length(x)` matrix
under mixed boundary conditions from `bc` using finite difference method specified by `method`. 
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   0.0   ⋅
 -1.0   1.0  0.0
   ⋅   -1.0  1.0

julia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅
  0.0  -1.0   1.0
   ⋅    0.0  -1.0

julia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 0.0   1.0    ⋅
 1.0  -2.0   1.0
  ⋅    1.0  -2.0
```
"""
function DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)
    T = eltype(x)
    d = diff(x)
    Δ_1 = d[1]
    Δ_M = d[end]
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ

    # get basis operator on interior nodes
    L = DifferentialOperator(x, (Reflecting(), Reflecting()), method)

    # apply boundary conditions
    L[1,1] -= typeof(method) <: BackwardFirstDifference ? ξ_lb : zero(T)
    L[1,1] += typeof(method) <: CentralSecondDifference ? (ξ_lb / Δ_1) : zero(T)
    L[end,end] -= typeof(method) <: ForwardFirstDifference ? ξ_ub : zero(T)
    L[end,end] -= typeof(method) <: CentralSecondDifference ? (ξ_ub / Δ_M) : zero(T) 

    return L
end

# Convenience calls
"""
    `L₁₋(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`
Returns a discretized first-order differential operator of `length(x)` by `length(x)` matrix
using backward difference under boundary conditions specified by `bc`.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper. 
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> L₁₋(x, (Reflecting(), Reflecting()))
3×3 Tridiagonal{Float64,Array{Float64,1}}:
  0.0   0.0   ⋅
 -1.0   1.0  0.0
   ⋅   -1.0  1.0
```
"""
L₁₋(x, bc) = DifferentialOperator(x, bc, BackwardFirstDifference())

"""
    `L₁₊(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`
Returns a discretized first-order differential operator of `length(x)` by `length(x)` matrix
using forward difference under boundary conditions specified by `bc`.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper. 
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> L₁₊(x, (Reflecting(), Reflecting()))
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0   ⋅
  0.0  -1.0  1.0
   ⋅    0.0  0.0
```
"""
L₁₊(x, bc) = DifferentialOperator(x, bc, ForwardFirstDifference())

"""
    `L₂(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`
Returns a discretized second-order differential operator of `length(x)` by `length(x)` matrix
using central difference under boundary conditions specified by `bc`.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper. 
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> L₂(x, (Reflecting(), Reflecting()))
3×3 Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅
  1.0  -2.0   1.0
   ⋅    1.0  -1.0
```
"""
L₂(x, bc) = DifferentialOperator(x, bc, CentralSecondDifference())

"""
    `L̄₁₋(x)`
Returns a discretized first-order differential operator of `length(x)` by `length(x) + 2` matrix
using backward difference under no boundary condition.

The first and last columns are applied to the ghost nodes just before `x[1]` and `x[end]` respectively.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> Array(L̄₁₋(x))
3×5 Array{Float64,2}:
 -1.0   1.0   0.0  0.0  0.0
  0.0  -1.0   1.0  0.0  0.0
  0.0   0.0  -1.0  1.0  0.0
```
"""
L̄₁₋(x) = ExtensionDifferentialOperator(x, BackwardFirstDifference())

"""
    `L̄₁₊(x)`
Returns a discretized first-order differential operator of `length(x)` by `length(x) + 2` matrix using 
forward difference under no boundary condition.

The first and last columns are applied to the ghost nodes just before `x[1]` and `x[end]` respectively.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> Array(L̄₁₊(x))
3×5 Array{Float64,2}:
 0.0  -1.0   1.0   0.0  0.0
 0.0   0.0  -1.0   1.0  0.0
 0.0   0.0   0.0  -1.0  1.0
```
"""
L̄₁₊(x) = ExtensionDifferentialOperator(x, ForwardFirstDifference())

"""
    `L̄₂(x)`
Returns a discretized second-order differential operator of `length(x)` by `length(x) + 2` matrix
using central difference under no boundary condition.

The first and last columns are applied to the ghost nodes just before `x[1]` and `x[end]` respectively.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> Array(L̄₂(x))
3×5 Array{Float64,2}:
 1.0  -2.0   1.0   0.0  0.0
 0.0   1.0  -2.0   1.0  0.0
 0.0   0.0   1.0  -2.0  1.0
```
"""
L̄₂(x)  = ExtensionDifferentialOperator(x, CentralSecondDifference())

"""
    `x̄(x)`
Returns an extended grid of length `length(x)+2` given grid `x`.

The first and last elements of the returned extended grid represent the ghost nodes
just before `x[1]` and `x[end]` respectively.
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> x̄(x)
5-element Array{Int64,1}:
 0
 1
 2
 3
 4

julia> x = [1.0; 1.5; 1.7]
3-element Array{Float64,1}:
 1.0
 1.5
 1.7

julia> x̄(x)
5-element Array{Float64,1}:
 0.5
 1.0
 1.5
 1.7
 1.9
```
"""
function x̄(x)
    d = diff(x) # dispatches based on AbstractArray or not
    x̄ = collect([x[1] - d[1]; x; x[end] + d[end]])
end

"""
    `diffusionoperators(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`
Returns a tuple of differential operators and extended grid `(L₁₋, L₁₊, L₂, x̄)`
with specified boundary conditions.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper. 
`x̄` is a `(M+2)` array that represents the extended grid whose 
first and last elements represent the ghost nodes just before `x[1]` and after `x[end]`.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> operators = diffusionoperators(x, (Reflecting(), Reflecting()));

julia> Array(operators.L₁₋)
3×3 Array{Float64,2}:
  0.0   0.0  0.0
 -1.0   1.0  0.0
  0.0  -1.0  1.0

julia> Array(operators.L₁₊)
3×3 Array{Float64,2}:
 -1.0   1.0  0.0
  0.0  -1.0  1.0
  0.0   0.0  0.0

julia> Array(operators.L₂)
3×3 Array{Float64,2}:
 -1.0   1.0   0.0
  1.0  -2.0   1.0
  0.0   1.0  -1.0

julia> operators.x̄
5-element Array{Int64,1}:
 0
 1
 2
 3
 4
```
"""
diffusionoperators(x, bc) = (L₁₋ = L₁₋(x, bc), L₁₊ = L₁₊(x, bc), L₂ = L₂(x, bc), x̄ = x̄(x))

"""
    `diffusionoperators(x)`
Returns a tuple of differential operators on extended grid and the extended grid itself `(L̄₁₋, L̄₁₊, L̄₂, x̄)`
without boundary conditions.

Given a grid `x` of length `M`, return diffusion operators for negative drift, positive drift,
and central differences of `M+2` by `M` matrices without boundary conditions. The first column is applied
for the ghost node just before `x[1]` and the last column for the one after `x[end]`.
`x̄` is a `(M+2)` array that represents the extended grid whose 
first and last elements represent the ghost nodes just before `x[1]` and after `x[end]`.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x = 1:3
1:3

julia> operators = diffusionoperators(x);

julia> Array(operators.L̄₁₋)
3×5 Array{Float64,2}:
 -1.0   1.0   0.0  0.0  0.0
  0.0  -1.0   1.0  0.0  0.0
  0.0   0.0  -1.0  1.0  0.0

julia> Array(operators.L̄₁₊)
3×5 Array{Float64,2}:
 0.0  -1.0   1.0   0.0  0.0
 0.0   0.0  -1.0   1.0  0.0
 0.0   0.0   0.0  -1.0  1.0

julia> Array(operators.L̄₂)
3×5 Array{Float64,2}:
 1.0  -2.0   1.0   0.0  0.0
 0.0   1.0  -2.0   1.0  0.0
 0.0   0.0   1.0  -2.0  1.0

julia> operators.x̄
5-element Array{Int64,1}:
 0
 1
 2
 3
 4
```
"""
diffusionoperators(x) = (L̄₁₋ = L̄₁₋(x), L̄₁₊ = L̄₁₊(x), L̄₂ = L̄₂(x), x̄ = x̄(x))
