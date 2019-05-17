function jointoperator_bc(operators, Q)
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
    d = reduce(vcat, [operator[band(0)] for operator in operators])
    dl = reduce(vcat, [[operator[band(-1)]; 0.] for operator in operators])[1:end-1]
    du = reduce(vcat, [[operator[band(1)]; 0.] for operator in operators])[1:end-1]
    Ls = BandedMatrix((-1 => dl, 0 => d, 1 => du), (M*N,M*N))
    Ls = BandedBlockBandedMatrix(Ls, (M*ones(Int64, N), M*ones(Int64, N)), (0,0), (1,1))

    # construct a kronecker product of Q times I_M
    Qs = BandedBlockBandedMatrix(Kron(Q, Eye(M)))

    return (Ls + Qs)
end
