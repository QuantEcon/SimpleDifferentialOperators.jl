"""
    L₀affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})

    Returns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`
to solve the system `L₀ * v(̄x) = b` on interior nodes.

    Returns a `length(interiornodes(x̄))`-length vector such that solving 
`L₀ * v(x̄) = b` on interior nodes `x` under the boundary condition given by `bc`
is identical with `L₀bc * v(x) = b + L₀affine(x̄, bc)` The first element of `bc` is 
applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> b = 1:3
1:3

julia> L₀affine(x̄, b, (Reflecting(), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> L₀affine(x̄, b, (NonhomogeneousAbsorbing(1.0), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> L₀affine(x̄, b, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
```
"""
function L₀affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    return zeros(M)
end

"""
    L₁₋affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})

    Returns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`
to solve the system `L₁₋ * v(̄x) = b` on interior nodes.

    Returns a `length(interiornodes(x̄))`-length vector such that solving 
`L₁₋ * v(x̄) = b` on interior nodes `x` under the boundary condition given by `bc`
is identical with `L₁₋bc * v(x) = b + L₁₋affine(x̄, bc)`. The first element of `bc` is 
applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> b = 1:3
1:3

julia> L₁₋affine(x̄, b, (Reflecting(), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> L₁₋affine(x̄, b, (NonhomogeneousAbsorbing(1.0), Reflecting()))
4-element Array{Float64,1}:
 1.0
 0.0
 0.0
 0.0

julia> L₁₋affine(x̄, b, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))
4-element Array{Float64,1}:
 1.0
 0.0
 0.0
 0.0
```
"""
function L₁₋affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    b_affine = zeros(M)

    # apply affine operator if needed
    Δ_1m = x̄[2] - x̄[1]
    b_affine[1] = (typeof(bc[1]) <: NonhomogeneousAbsorbing) ? bc[1].S / Δ_1m : b_affine[1]

    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        b_affine[1:min(M, bc[1].loc)] = -b[1:min(M, bc[1].loc)]
    end
    
    return b_affine
end
"""
    L₁₊affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})

    Returns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`
to solve the system `L₁₊ * v(̄x) = b` on interior nodes.

    Returns a `length(interiornodes(x̄))`-length vector such that solving 
`L₁₊ * v(x̄) = b` on interior nodes `x` under the boundary condition given by `bc`
is identical with `L₁₊bc * v(x) = b + L₁₊affine(x̄, bc)` The first element of `bc` is 
applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> b = 1:3
1:3

julia> L₁₊affine(x̄, b, (Reflecting(), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> L₁₊affine(x̄, b, (NonhomogeneousAbsorbing(1.0), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> L₁₊affine(x̄, b, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))
4-element Array{Float64,1}:
  0.0
  0.0
  0.0
 -2.0
```
"""
function L₁₊affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    b_affine = zeros(M)

    # apply affine operator if needed
    Δ_Mp = x̄[end] - x̄[end-1]
    b_affine[end] = (typeof(bc[2]) <: NonhomogeneousAbsorbing) ? -bc[2].S / Δ_Mp : b_affine[end]

    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        b_affine[1:min(M, bc[1].loc)] = -b[1:min(M, bc[1].loc)]
    end
    
    return b_affine
end

"""
    L₂affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})

    Returns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`
to solve the system `L₂ * v(̄x) = b` on interior nodes.

    Returns a `length(interiornodes(x̄))`-length vector such that solving 
`L₂ * v(x̄) = b` on interior nodes `x` under the boundary condition given by `bc`
is identical with `L₂bc * v(x) = b + L₂affine(x̄, bc)` The first element of `bc` is 
applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> b = 1:3
1:3

julia> L₂affine(x̄, b, (Reflecting(), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> L₂affine(x̄, b, (NonhomogeneousAbsorbing(1.0), Reflecting()))
4-element Array{Float64,1}:
-1.0
 0.0
 0.0
 0.0

julia> L₂affine(x̄, b, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))
4-element Array{Float64,1}:
-1.0
 0.0
 0.0
-2.0
```
"""
function L₂affine(x̄, b, bc::Tuple{BoundaryCondition, BoundaryCondition})
    x = interiornodes(x̄)
    M = length(x)
    b_affine = zeros(M)

    # apply affine operator if needed
    Δ_1p = x̄[3] - x̄[2]
    Δ_1m = x̄[2] - x̄[1]
    Δ_Mp = x̄[end] - x̄[end-1]
    Δ_Mm = x̄[end-1] - x̄[end-2]
    b_affine[1] = (typeof(bc[1]) <: NonhomogeneousAbsorbing) ? -2*bc[1].S/(Δ_1m*(Δ_1m+Δ_1p)) : b_affine[1]
    b_affine[end] = (typeof(bc[2]) <: NonhomogeneousAbsorbing) ? -2*bc[2].S/(Δ_Mp*(Δ_Mm+Δ_Mp)) : b_affine[end]

    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        b_affine[1:min(M, bc[1].loc)] = -b[1:min(M, bc[1].loc)]
    end

    return b_affine
end

"""
    Laffine(L, bc::Tuple{BoundaryCondition, BoundaryCondition})

    Returns a vector for affine operator given an extension operator `L` and boundary conditions `bc`
to solve the system `L * v(̄x) = b` on interior nodes.

    Returns a `size(L)[1]`-length vector such that solving 
`L * v(x̄) = b` on interior nodes `x` under the boundary condition given by `bc`
is identical with `Lbc * v(x) = b + Laffine(L, bc)` The first element of `bc` is 
applied to the lower bound, and second element of `bc` to the upper.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> b = 1:3
1:3

julia> L = L₂(x̄)
4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 12 stored entries:
  [1, 1]  =  1.0
  [1, 2]  =  -2.0
  [2, 2]  =  1.0
  [1, 3]  =  1.0
  [2, 3]  =  -2.0
  [3, 3]  =  1.0
  [2, 4]  =  1.0
  [3, 4]  =  -2.0
  [4, 4]  =  1.0
  [3, 5]  =  1.0
  [4, 5]  =  -2.0
  [4, 6]  =  1.0

julia> Laffine(L, b, (Reflecting(), Reflecting()))
4-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0

julia> Laffine(L, b, (NonhomogeneousAbsorbing(1.0), Reflecting()))
4-element Array{Float64,1}:
-1.0
 0.0
 0.0
 0.0

julia> Laffine(L, b, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))
4-element Array{Float64,1}:
-1.0
 0.0
 0.0
-2.0
```
"""
# take M by (M+2) extended operator 
function Laffine(L, b, bc::Tuple{BoundaryCondition, BoundaryCondition})
    M = size(L)[1]
    b_affine = zeros(M)
    
    # apply affine operator if needed
    b_affine[1] = (typeof(bc[1]) <: NonhomogeneousAbsorbing) ? -sum(L[:,1])*bc[1].S : b_affine[1]
    b_affine[end] = (typeof(bc[2]) <: NonhomogeneousAbsorbing) ? -sum(L[:,end])*bc[2].S : b_affine[end]
    
    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        b_affine[1:min(M, bc[1].loc)] = -b[1:min(M, bc[1].loc)]
    end
    
    return b_affine
end