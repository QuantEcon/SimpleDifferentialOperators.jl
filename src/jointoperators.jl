"""
    jointoperator_bc(operators, Q::Array)

Returns a discretized operator that solves systems of differential equations defined by
`operators` with transitions by `Q` where `operators` is an N-length collection of 
discretized operators with boundary conditions applied and `Q` is N by N intensity matrix.

# Examples
```jldoctest; setup = :(using SimpleDifferentialOperators)
julia> x̄ = 0:5
0:5

julia> L1bc = L₁₋bc(x̄, (Reflecting(), Reflecting()))
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
  0.0   0.0    ⋅    ⋅
 -1.0   1.0   0.0   ⋅
   ⋅   -1.0   1.0  0.0
   ⋅     ⋅   -1.0  1.0

julia> L2bc = L₂bc(x̄, (Reflecting(), Reflecting()))
4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:
 -1.0   1.0    ⋅     ⋅
  1.0  -2.0   1.0    ⋅
   ⋅    1.0  -2.0   1.0
   ⋅     ⋅    1.0  -1.0

julia> Q = [-0.5 0.5; 0.3 -0.3]
2×2 Array{Float64,2}:
 -0.5   0.5
  0.3  -0.3

julia> jointoperator_bc((L1bc, L2bc), Q)
2×2-blocked 8×8 BlockBandedMatrices.BandedBlockBandedMatrix{Float64,BlockArrays.PseudoBlockArray{Float64,2,Array{Float64,2},BlockArrays.BlockSizes{2,Tuple{Array{Int64,1},Array{Int64,1}}}}}:
 -0.5   0.0    ⋅    ⋅   │   0.5   0.0    ⋅     ⋅
 -1.0   0.5   0.0   ⋅   │   0.0   0.5   0.0    ⋅
  ⋅   -1.0   0.5  0.0  │    ⋅    0.0   0.5   0.0
  ⋅     ⋅   -1.0  0.5  │    ⋅     ⋅    0.0   0.5
───────────────────────┼────────────────────────
  0.3   0.0    ⋅    ⋅   │  -1.3   1.0    ⋅     ⋅
  0.0   0.3   0.0   ⋅   │   1.0  -2.3   1.0    ⋅
  ⋅    0.0   0.3  0.0  │    ⋅    1.0  -2.3   1.0
  ⋅     ⋅    0.0  0.3  │    ⋅     ⋅    1.0  -1.3
```
"""
function jointoperator_bc(operators, Q::Array)
    M = size(operators[1],1)
    N = length(operators)

    # check if all operators are square
    @assert all(operator->(size(operator,1) == size(operator,2)), 
                operators)

    # check if all operators have same size
    @assert all(operator->(size(operator) == size(operators[1])), 
                operators)

    # check if the size of transition matrix is 
    # same as the number of operators 
    @assert size(Q,1) == size(Q,2) == N

    # extract operators and append them to form a diagonal block tridiagonal 
    Ls = blockdiag(sparse.(operators)...)
    Ls = BandedBlockBandedMatrix(Ls, (M*ones(Int64, N), M*ones(Int64, N)), (0,0), (1,1))

    # construct a kronecker product of Q times I_M
    Qs = BandedBlockBandedMatrix(Kron(Q, Eye(M)))

    return (Ls + Qs)
end
