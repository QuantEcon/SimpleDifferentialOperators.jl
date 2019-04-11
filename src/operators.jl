"""
    DifferentialOperator(x̄, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)

Returns a discretized differential operator of `length(x̄)` by `length(x̄)` matrix
under reflecting boundary conditions from `bc` using finite difference method specified by `method`.

# Examples

```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), BackwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
  0.0   0.0    ⋅    ⋅
 -1.0   1.0   0.0   ⋅
   ⋅   -1.0   1.0  0.0
   ⋅     ⋅   -1.0  1.0

julia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), ForwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅    ⋅
  0.0  -1.0   1.0   ⋅
   ⋅    0.0  -1.0  1.0
   ⋅     ⋅    0.0  0.0

julia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), CentralSecondDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅     ⋅
  1.0  -2.0   1.0    ⋅
   ⋅    1.0  -2.0   1.0
   ⋅     ⋅    1.0  -1.0
```
"""
function DifferentialOperator(x̄, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)
    T = eltype(x̄)
    Δ_1m = x̄[2] - x̄[1]
    Δ_1p = x̄[3] - x̄[2]
    Δ_Mm = x̄[end-1] - x̄[end-2]
    Δ_Mp = x̄[end] - x̄[end-1]

    # get basis operator on interior nodes
    L = get_basis_operator(x̄, method)
    
    Ξ_1 = -(1/(Δ_1m*Δ_1p)-1/((Δ_1p+Δ_1m)*(Δ_1m)))
    Ξ_M = -(1/(Δ_Mm*Δ_Mp)-1/((Δ_Mp+Δ_Mm)*(Δ_Mp))) 
   
    # apply boundary conditions
    L[1,1] = typeof(method) <: BackwardFirstDifference ? zero(T) : L[1,1]
    L[1,1] = typeof(method) <: CentralSecondDifference ? 2*Ξ_1 : L[1,1]
    L[end,end] = typeof(method) <: ForwardFirstDifference ? zero(T) : L[end,end]
    L[end,end] = typeof(method) <: CentralSecondDifference ? 2*Ξ_M : L[end,end]

    return L
end

"""
    DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)

Returns a discretized differential operator of `length(x̄)` by `length(x̄)` matrix
under mixed boundary conditions from `bc` using finite difference method specified by `method`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   0.0    ⋅    ⋅
 -1.0   1.0   0.0   ⋅
   ⋅   -1.0   1.0  0.0
   ⋅     ⋅   -1.0  1.0

julia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅     ⋅
  0.0  -1.0   1.0    ⋅
   ⋅    0.0  -1.0   1.0
   ⋅     ⋅    0.0  -1.0

julia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 0.0   1.0    ⋅     ⋅
 1.0  -2.0   1.0    ⋅
  ⋅    1.0  -2.0   1.0
  ⋅     ⋅    1.0  -2.0
```
"""
function DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)
    T = eltype(x̄)
    d = diff(x̄)
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ
    Δ_1p = x̄[3] - x̄[2]
    Δ_1m = x̄[2] - x̄[1]
    Δ_Mp = x̄[end] - x̄[end-1]
    Δ_Mm = x̄[end-1] - x̄[end-2]

    # get extended operator with reflecting barrier conditions first
    L = get_basis_operator(x̄, method)
 
    Ξ_1 = -(1/(Δ_1m*Δ_1p)+1/((-1+ξ_lb*Δ_1m)*(Δ_1p+Δ_1m)*(Δ_1m)))
    Ξ_M = -(1/(Δ_Mm*Δ_Mp)-1/((1+ξ_ub*Δ_Mp)*(Δ_Mp+Δ_Mm)*(Δ_Mp)))
   
    # apply boundary conditions
    L[1,1] = typeof(method) <: BackwardFirstDifference ? (1+1/(-1+ξ_lb*Δ_1m))/Δ_1m : L[1,1]
    L[1,1] = typeof(method) <: CentralSecondDifference ? 2*Ξ_1 : L[1,1]
    L[end,end] = typeof(method) <: ForwardFirstDifference ? (-1+1/(1+ξ_ub*Δ_Mp))/Δ_Mp : L[end,end]
    L[end,end] = typeof(method) <: CentralSecondDifference ? 2*Ξ_M : L[end,end]

    return L
end

# Convenience calls
"""
    L₁₋bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})

Returns a discretized first-order differential operator of `length(x̄)` by `length(x̄)` matrix
using backward difference under boundary conditions specified by `bc`.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> L₁₋bc(x̄, (Reflecting(), Reflecting()))
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
  0.0   0.0    ⋅    ⋅
 -1.0   1.0   0.0   ⋅
   ⋅   -1.0   1.0  0.0
   ⋅     ⋅   -1.0  1.0
```
"""
L₁₋bc(x̄, bc) = DifferentialOperator(x̄, bc, BackwardFirstDifference())

"""
    L₁₊bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})

Returns a discretized first-order differential operator of `length(x̄)` by `length(x̄)` matrix
using forward difference under boundary conditions specified by `bc`.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> L₁₊bc(x̄, (Reflecting(), Reflecting()))
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅    ⋅
  0.0  -1.0   1.0   ⋅
   ⋅    0.0  -1.0  1.0
   ⋅     ⋅    0.0  0.0
```
"""
L₁₊bc(x̄, bc) = DifferentialOperator(x̄, bc, ForwardFirstDifference())

"""
    L₂bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})

Returns a discretized second-order differential operator of `length(x̄)` by `length(x̄)` matrix
using central difference under boundary conditions specified by `bc`.

The first element of `bc` is applied to the lower bound, and second element of `bc` to the upper.
# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> L₂bc(x̄, (Reflecting(), Reflecting()))
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅     ⋅
  1.0  -2.0   1.0    ⋅
   ⋅    1.0  -2.0   1.0
   ⋅     ⋅    1.0  -1.0
```
"""
L₂bc(x̄, bc) = DifferentialOperator(x̄, bc, CentralSecondDifference())

"""
    L₁₋(x̄)

Returns a discretized first-order differential operator of `length(x̄)` by `length(x̄) + 2` matrix
using backward difference under no boundary condition.

The first and last columns are applied to the ghost nodes just before `x̄[1]` and `x̄[end]` respectively.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 1:3
1:3

julia> Array(L₁₋(x̄))
1×3 Array{Float64,2}:
 -1.0  1.0  0.0
```
"""
L₁₋(x̄) = ExtensionDifferentialOperator(x̄, BackwardFirstDifference())

"""
    L₁₊(x̄)

Returns a discretized first-order differential operator of `length(x̄)` by `length(x̄) + 2` matrix using
forward difference under no boundary condition.

The first and last columns are applied to the ghost nodes just before `x̄[1]` and `x̄[end]` respectively.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> Array(L₁₊(x̄))
4×6 Array{Float64,2}:
 0.0  -1.0   1.0   0.0   0.0  0.0
 0.0   0.0  -1.0   1.0   0.0  0.0
 0.0   0.0   0.0  -1.0   1.0  0.0
 0.0   0.0   0.0   0.0  -1.0  1.0
```
"""
L₁₊(x̄) = ExtensionDifferentialOperator(x̄, ForwardFirstDifference())

"""
    L₂(x̄)

Returns a discretized second-order differential operator of `length(x̄)` by `length(x̄) + 2` matrix
using central difference under no boundary condition.

The first and last columns are applied to the ghost nodes just before `x̄[1]` and `x̄[end]` respectively.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5 

julia> Array(L₂(x̄))
4×6 Array{Float64,2}:
 1.0  -2.0   1.0   0.0   0.0  0.0
 0.0   1.0  -2.0   1.0   0.0  0.0
 0.0   0.0   1.0  -2.0   1.0  0.0
 0.0   0.0   0.0   1.0  -2.0  1.0
```
"""
L₂(x̄)  = ExtensionDifferentialOperator(x̄, CentralSecondDifference())

"""
    interiornodes(x̄)

Returns an interior grid of length `length(x̄)-2` given extended grid `x̄`.
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> interiornodes(x̄)
1:4

julia> x̄ = [1.0; 1.5; 1.7]
3-element Array{Float64,1}:
 1.0
 1.5
 1.7

julia> interiornodes(x̄)
1-element Array{Float64,1}:
 1.5
```
"""
interiornodes(x̄) = x̄[2:end-1]

"""
    interiornodes(x̄, bc)

Returns an interior grid corresponding to the boundary condition `bc` given extended grid `x̄`.
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> interiornodes(x̄, (Reflecting(), Reflecting()))
1:4

julia> x̄ = [1.0; 1.5; 1.7]
3-element Array{Float64,1}:
 1.0
 1.5
 1.7

julia> interiornodes(x̄, (Mixed(1.0), Mixed(1.0)))
1-element Array{Float64,1}:
 1.5
```
"""
interiornodes(x̄, bc) = interiornodes(x̄)