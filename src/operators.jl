"""
    DifferentialOperator(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition}, method::DiscretizationMethod)

Returns a discretized differential operator of 
`length(interiornodes(x̄))` by `length(interiornodes(x̄))` matrix
under mixed boundary conditions from `bc` using a discretization method specified by `method`.

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

julia> x̄ = 0:5
0:5

julia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), BackwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 Inf     0.0    ⋅    ⋅
  -1.0   1.0   0.0   ⋅
    ⋅   -1.0   1.0  0.0
    ⋅     ⋅   -1.0  1.0

julia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), ForwardFirstDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅     ⋅
  0.0  -1.0   1.0    ⋅
   ⋅    0.0  -1.0   1.0
   ⋅     ⋅    0.0  -0.5

julia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), CentralSecondDifference())
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -Inf     1.0    ⋅     ⋅
    1.0  -2.0   1.0    ⋅
     ⋅    1.0  -2.0   1.0
     ⋅     ⋅    1.0  -1.5
```
"""
function DifferentialOperator(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition}, method::BackwardFirstDifference)
    # reflecting bcs are special cases of mixed bcs with ξ = 0
    if (typeof(bc[1]) <: Reflecting)
        return DifferentialOperator(x̄, (Mixed(ξ = 0), bc[2]), method)
    end
    if (typeof(bc[2]) <: Reflecting)
        return DifferentialOperator(x̄, (bc[1], Mixed(ξ = 0)), method)
    end

    # setup for operator
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₋⁻¹ = 1 ./ d[1:end-1] # (1 ./ Δ₋)
    T = eltype(Δ₋⁻¹)

    # construct the operator
    L = Tridiagonal(-Δ₋⁻¹[2:M], Δ₋⁻¹, zeros(T, M-1)) 

    # setup for boundary conditions
    Δ_1m = x̄[2] - x̄[1]

    # apply boundary conditions
    # (under homogeneous absorbing bc on lb, the first column in invariant)
    if (typeof(bc[1]) <: Mixed) 
        ξ_lb = bc[1].ξ
        L[1,1] += (bc[1].direction == :backward) ? (-1/Δ_1m - ξ_lb) : 1/(-1+ξ_lb*Δ_1m)/Δ_1m
    end
    # apply absorbing boundary condition on lb further if it's applied withint boundary
    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        L[:,1:(min(M, bc[1].loc))] .= zero(T)
    end

    return L
end

function DifferentialOperator(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition}, method::ForwardFirstDifference)
    # reflecting bcs are special cases of mixed bcs with ξ = 0
    if (typeof(bc[1]) <: Reflecting)
        return DifferentialOperator(x̄, (Mixed(ξ = 0), bc[2]), method)
    end
    if (typeof(bc[2]) <: Reflecting)
        return DifferentialOperator(x̄, (bc[1], Mixed(ξ = 0)), method)
    end

    # setup for operator
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₊⁻¹ = 1 ./ d[2:end] # (1 ./ Δ₊), extracting elements on the interior
    T = eltype(Δ₊⁻¹)

    # construct the operator
    L = Tridiagonal(zeros(T, M-1), -Δ₊⁻¹, Δ₊⁻¹[1:M-1]) 
    
    # setup for boundary conditions
    Δ_Mp = x̄[end] - x̄[end-1]

    # apply boundary conditions
    # (under homogeneous absorbing bc on ub, the last column in invariant)
    if (typeof(bc[2]) <: Mixed) 
        ξ_ub = bc[2].ξ
        L[end,end] += (bc[2].direction == :forward) ? (1/Δ_Mp - ξ_ub) : 1/(1+ξ_ub*Δ_Mp)/Δ_Mp
    end
    
    # apply absorbing boundary condition on lb further if it's applied withint boundary
    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        L[:,1:(min(M, bc[1].loc))] .= zero(T)
    end

    return L
end

function DifferentialOperator(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition}, method::CentralSecondDifference)
    # reflecting bcs are special cases of mixed bcs with ξ = 0
    if (typeof(bc[1]) <: Reflecting)
        return DifferentialOperator(x̄, (Mixed(ξ = 0), bc[2]), method)
    end
    if (typeof(bc[2]) <: Reflecting)
        return DifferentialOperator(x̄, (bc[1], Mixed(ξ = 0)), method)
    end

    # setup for operators
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₋⁻¹ = 1 ./ d[1:end-1] # 1 ./ Δ₋
    Δ₊⁻¹ = 1 ./ d[2:end] # 1 ./ Δ₊
    Δ⁻¹ = 1 ./ (d[1:end-1] + d[2:end]) # 1 ./ (Δ₋ + Δ₊)
    T = eltype(Δ⁻¹)

    # construct the operator
    L = 2*Tridiagonal((Δ⁻¹.*Δ₋⁻¹)[2:M], -Δ₋⁻¹ .* Δ₊⁻¹, (Δ⁻¹.*Δ₊⁻¹)[1:M-1])

    # setup for boundary conditions
    Δ_1p = x̄[3] - x̄[2]
    Δ_1m = x̄[2] - x̄[1]
    Δ_Mp = x̄[end] - x̄[end-1]
    Δ_Mm = x̄[end-1] - x̄[end-2]

    # apply boundary conditions
    # (under homogeneous absorbing bc on lb, the first column in invariant)
    if (typeof(bc[1]) <: Mixed) 
        ξ_lb = bc[1].ξ
        Ξ_1p = L[1,1] - 2/((-1+ξ_lb*Δ_1m)*(Δ_1p+Δ_1m)*(Δ_1m))
        Ξ_1m = 2*(-1/(Δ_1p*Δ_1m) + (1+ξ_lb*Δ_1m)/(Δ_1p+Δ_1m)/Δ_1m)
        L[1,1] = (bc[1].direction == :backward) ? Ξ_1m : Ξ_1p
    end
    # (under homogeneous absorbing bc on ub, the last column in invariant)
    if (typeof(bc[2]) <: Mixed) 
        ξ_ub = bc[2].ξ
        Ξ_Mm = L[end,end] + 2/((1+ξ_ub*Δ_Mp)*(Δ_Mp+Δ_Mm)*(Δ_Mp))
        Ξ_Mp = 2*(-1/(Δ_Mp*Δ_Mm) - (-1+ξ_ub*Δ_Mp)/(Δ_Mp+Δ_Mm)/Δ_Mp)    
        L[end,end] = (bc[2].direction == :forward) ? Ξ_Mp : Ξ_Mm
    end

    # apply absorbing boundary condition on lb further if it's applied withint boundary
    if (typeof(bc[1]) <: Absorbing && bc[1].loc > 1)
        L[:,1:(min(M, bc[1].loc))] .= zero(T)
    end
    
    return L
end

function JumpProcessOperator(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition}, method::JumpProcess)
    T = eltype(x̄)
    M = length(x̄) - 2
    jumps = method.jumps
    destinations = Array(1:M) + jumps

    # under reflecting or mixed boundary conditions, 
    # we assume that there is no destination on the boundary nodes 
    if ((typeof(bc[1]) <: Reflecting || typeof(bc[1]) <: Mixed) && any(destinations .<= 0))
        throw(ErrorException("Under reflecting/mixed boundary conditions on lower boundary, the lower boundary node cannot be destination."))
    end
    if ((typeof(bc[2]) <: Reflecting || typeof(bc[2]) <: Mixed) && any(destinations .>= (M+1))) 
        throw(ErrorException("Under reflecting/mixed boundary conditions on upper boundary, the upper boundary node cannot be destination."))
    end

    # construct L
    L = BandedMatrix((0=>-ones(M), 1=>zeros(M-1)), (M,M), 
                    (max(abs(minimum(jumps)), 0), max(abs(maximum(jumps)), 1)))

    # apply destination for each row
    for i in 1:length(jumps)
        # add destination if it is not on/out of the boundary
        if !(destinations[i] <= 0 || destinations[i] >= M+1)
            L[i,destinations[i]] += 1
        end
    end

    return L
end


# Convenience calls
"""
    L₁₋bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})

Returns a discretized first-order differential operator of 
`length(interiornodes(x̄))` by `length(interiornodes(x̄))` matrix
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

Returns a discretized first-order differential operator of 
`length(interiornodes(x̄))` by `length(interiornodes(x̄))` matrix
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

Returns a discretized second-order differential operator of 
`length(interiornodes(x̄))` by `length(interiornodes(x̄))` matrix
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
    Lₙbc(x̄, method, (Absorbing(), Absorbing()))

Returns a discretized jump process operator of
`length(interiornodes(x̄))` by `length(interiornodes(x̄))` matrix
specified by `method` under boundary conditions specified by `bc`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5 

julia> Lₙbc(x̄, (Absorbing(), Absorbing()), JumpProcess(x̄, -1.0))
4×4 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0   0.0    ⋅     ⋅
 1.0  -1.0   0.0    ⋅
  ⋅    1.0  -1.0   0.0
  ⋅     ⋅    1.0  -1.0
```
"""
Lₙbc(x̄, bc, method)  = JumpProcessOperator(x̄, bc, method)

"""
    L₁₋(x̄)

Returns a discretized first-order differential operator of 
`length(interiornodes(x̄))` by `length(x̄)` matrix
using backward difference under no boundary condition.

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

Returns a discretized first-order differential operator of
`length(interiornodes(x̄))` by `length(x̄)` matrix
using forward difference under no boundary condition.

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

Returns a discretized second-order differential operator of
`length(interiornodes(x̄))` by `length(x̄)` matrix
using central difference under no boundary condition.

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
    Lₙ(x̄, method)

Returns a discretized jump process operator of
`length(interiornodes(x̄))` by `length(x̄)` matrix
specified by `method`

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5 

julia> Lₙ(x̄, JumpProcess(x̄, -1.0))
4×6 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:
 0.0  0.0    ⋅     ⋅     ⋅    ⋅
  ⋅   1.0  -1.0    ⋅     ⋅    ⋅
  ⋅    ⋅    1.0  -1.0    ⋅    ⋅
  ⋅    ⋅     ⋅    1.0  -1.0   ⋅
```
"""
Lₙ(x̄, method)  = ExtensionDifferentialOperator(x̄, method)

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