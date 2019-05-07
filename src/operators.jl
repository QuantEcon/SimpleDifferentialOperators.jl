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
function DifferentialOperator(x̄, bc::Tuple{Reflecting, Reflecting}, method)
    # since reflecting boundary conditions are special cases of
    # mixed boundary conditions with ξ_lb = ξ_ub = 0
    DifferentialOperator(x̄, (Mixed(ξ = 0), Mixed(ξ = 0)), method)
end

"""
    DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)

Returns a discretized differential operator of `length(x̄)` by `length(x̄)` matrix
under mixed boundary conditions from `bc` using finite difference method specified by `method`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), BackwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   0.0    ⋅    ⋅
 -1.0   1.0   0.0   ⋅
   ⋅   -1.0   1.0  0.0
   ⋅     ⋅   -1.0  1.0

julia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), ForwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅     ⋅
  0.0  -1.0   1.0    ⋅
   ⋅    0.0  -1.0   1.0
   ⋅     ⋅    0.0  -1.0

julia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), CentralSecondDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 0.0   1.0    ⋅     ⋅
 1.0  -2.0   1.0    ⋅
  ⋅    1.0  -2.0   1.0
  ⋅     ⋅    1.0  -2.0
```
"""
function DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::BackwardFirstDifference)
    # setup for operator
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₋⁻¹ = 1 ./ d[1:end-1] # (1 ./ Δ₋)
    T = eltype(Δ₋⁻¹)

    # construct the operator
    L = Tridiagonal(-Δ₋⁻¹[2:M], Δ₋⁻¹, zeros(T, M-1)) 

    # setup for boundary conditions
    ξ_lb = bc[1].ξ
    Δ_1m = x̄[2] - x̄[1]

    # apply boundary conditions
    L[1,1] += (bc[1].direction == :backward) ? (-1/Δ_1m - ξ_lb*Δ_1m) : 1/(-1+ξ_lb*Δ_1m)/Δ_1m

    return L
end

function DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::ForwardFirstDifference)
    # setup for operator
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₊⁻¹ = 1 ./ d[2:end] # (1 ./ Δ₊), extracting elements on the interior
    T = eltype(Δ₊⁻¹)

    # construct the operator
    L = Tridiagonal(zeros(T, M-1), -Δ₊⁻¹, Δ₊⁻¹[1:M-1]) 
    
    # setup for boundary conditions
    ξ_ub = bc[2].ξ
    Δ_Mp = x̄[end] - x̄[end-1]

    # apply boundary conditions
    L[end,end] += (bc[2].direction == :forward) ? (1/Δ_Mp - ξ_ub*Δ_Mp) : 1/(1+ξ_ub*Δ_Mp)/Δ_Mp

    return L
end

function DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::CentralSecondDifference)
    # setup for operators
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₋⁻¹ = 1 ./ d[1:end-1] # 1 ./ Δ₋
    Δ₊⁻¹ = 1 ./ d[2:end] # 1 ./ Δ₊
    Δ⁻¹ = 1 ./ (d[1:end-1] + d[2:end]) # 1 ./ (Δ₋ + Δ₊)

    # construct the operator
    L = 2*Tridiagonal((Δ⁻¹.*Δ₋⁻¹)[2:M], -Δ₋⁻¹ .* Δ₊⁻¹, (Δ⁻¹.*Δ₊⁻¹)[1:M-1])

    # setup for boundary conditions
    ξ_lb = bc[1].ξ
    ξ_ub = bc[2].ξ
    Δ_1p = x̄[3] - x̄[2]
    Δ_1m = x̄[2] - x̄[1]
    Δ_Mp = x̄[end] - x̄[end-1]
    Δ_Mm = x̄[end-1] - x̄[end-2]
    
    # apply boundary conditions
    Ξ_1p = L[1,1] - 2/((-1+ξ_lb*Δ_1m)*(Δ_1p+Δ_1m)*(Δ_1m))
    Ξ_Mm = L[end,end] + 2/((1+ξ_ub*Δ_Mp)*(Δ_Mp+Δ_Mm)*(Δ_Mp))
    Ξ_1m = 2*(-1/(Δ_1p*Δ_1m) + (1+ξ_lb*Δ_1m)/(Δ_1p+Δ_1m)/Δ_1m)
    Ξ_Mp = 2*(-1/(Δ_Mp*Δ_Mm) - (-1+ξ_ub*Δ_Mp)/(Δ_Mp+Δ_Mm)/Δ_Mp)

    L[1,1] = (bc[1].direction == :backward) ? Ξ_1m : Ξ_1p
    L[end,end] = (bc[2].direction == :forward) ? Ξ_Mp : Ξ_Mm
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

julia> interiornodes(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)))
1-element Array{Float64,1}:
 1.5
```
"""
interiornodes(x̄, bc) = interiornodes(x̄)