"""
    ExtensionDifferentialOperator(x̄, method::DifferenceMethod)

Returns a discretized differential operator of `length(x̄)` by `length(x̄) + 2` matrix
whose first and last columns are applied to the ghost nodes just before `x̄[1]` and `x̄[end]` respectively
under no boundary condition using finite difference method specified by `method`.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> ExtensionDifferentialOperator(x̄, BackwardFirstDifference())
4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 11 stored entries:
  [1, 1]  =  -1.0
  [1, 2]  =  1.0
  [2, 2]  =  -1.0
  [1, 3]  =  0.0
  [2, 3]  =  1.0
  [3, 3]  =  -1.0
  [2, 4]  =  0.0
  [3, 4]  =  1.0
  [4, 4]  =  -1.0
  [3, 5]  =  0.0
  [4, 5]  =  1.0

julia> ExtensionDifferentialOperator(x̄, ForwardFirstDifference())
4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 11 stored entries:
  [1, 2]  =  -1.0
  [2, 2]  =  0.0
  [1, 3]  =  1.0
  [2, 3]  =  -1.0
  [3, 3]  =  0.0
  [2, 4]  =  1.0
  [3, 4]  =  -1.0
  [4, 4]  =  0.0
  [3, 5]  =  1.0
  [4, 5]  =  -1.0
  [4, 6]  =  1.0

julia> ExtensionDifferentialOperator(x̄, CentralSecondDifference())
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
```
"""
function ExtensionDifferentialOperator(x̄::AbstractRange, method::BackwardFirstDifference)
    T = eltype(x̄)
    M = length(x̄) - 2
    Δ = step(x̄)
    return spdiagm(0 => -ones(T, M), 1 => ones(T, M), 2 => zeros(T, M))[1:M,:] / Δ
end

function ExtensionDifferentialOperator(x̄::AbstractRange, method::ForwardFirstDifference)
    T = eltype(x̄)
    M = length(x̄) - 2
    Δ = step(x̄)
    return spdiagm(0 => zeros(T, M), 1 => -ones(T, M), 2 => ones(T, M))[1:M,:] / Δ
end

function ExtensionDifferentialOperator(x̄::AbstractRange, method::CentralSecondDifference)
    T = eltype(x̄)
    M = length(x̄) - 2
    Δ = step(x̄)
    return spdiagm(0 => ones(T, M), 1 => -2*ones(T, M), 2 => ones(T, M))[1:M,:] / Δ^2
end

function ExtensionDifferentialOperator(x̄::AbstractArray, method::BackwardFirstDifference)
    T = eltype(x̄)
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₋⁻¹ = 1 ./ d[1:end-1]
    return spdiagm(0 => -Δ₋⁻¹, 1 => Δ₋⁻¹, 2 => zeros(T, M))[1:M,:]
end

function ExtensionDifferentialOperator(x̄::AbstractArray, method::ForwardFirstDifference)
    T = eltype(x̄)
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₊⁻¹ = 1 ./ d[2:end]
    return spdiagm(0 => zeros(T, M), 1 => -Δ₊⁻¹, 2 => Δ₊⁻¹)[1:M,:]
end

function ExtensionDifferentialOperator(x̄::AbstractArray, method::CentralSecondDifference)
    T = eltype(x̄)
    M = length(x̄) - 2
    d = diff(x̄)
    Δ₋⁻¹ = 1 ./ d[1:end-1] # 1 ./ Δ₋
    Δ₊⁻¹ = 1 ./ d[2:end] # 1 ./ Δ₊
    Δ⁻¹ = 1 ./ (d[1:end-1] + d[2:end]) # 1 ./ (Δ₋ + Δ₊)
    return 2*spdiagm(0 => Δ⁻¹.*Δ₋⁻¹, 1 => -Δ₋⁻¹ .* Δ₊⁻¹, 2 => Δ⁻¹.*Δ₊⁻¹)[1:M,:]
end

function ExtensionDifferentialOperator(x̄::AbstractArray, method::JumpProcess)
    T = eltype(x̄)
    M = length(x̄) - 2
    jumps = method.jumps

    L̄ = BandedMatrix((0=>zeros(M), 1=>-ones(M)), (M,M+2), 
                    (max(abs(minimum(jumps))-1, 0), max(abs(maximum(jumps)+1), 1)))
    for i in 1:length(jumps)
        L̄[i,(i+1+jumps[i])] += 1
    end

    return (L̄)
end