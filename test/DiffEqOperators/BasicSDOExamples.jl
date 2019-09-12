
using LinearAlgebra, DiffEqOperators, SimpleDifferentialOperators


# Below is a basic example from SDO:

# setup
π(x) = x^2
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of grid (interior points)

x̄ = range(0.0, 1.0, length = (M+2))
x = interiornodes(x̄) # i.e., x̄[2:end-1]

# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to reflecting barriers at 0 and 1
bc = (Reflecting(), Reflecting())
Lₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)
L_bc = I * ρ - Lₓ

# solve the value function
v = L_bc \ π.(x)


# The DiffEq approach (as of the current implementation)

# setup
π(x) = x^2
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of grid (interior points)

x̄ = range(0.0, 1.0, length = (M+2))
x = interiornodes(x̄) # i.e., x̄[2:end-1]
dx = x[2] - x[1]

# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to reflecting barriers at 0 and 1

L1 = UpwindDifference(1,1,dx,M,t->1.0)
L2 = CenteredDifference(2,2,dx,M)
Q = Neumann0BC(dx, 1)

# Here Array(A::GhostDerivativeOperator) will return a tuple of the linear part
# and the affine part of the operator A, hence we index Array(µ*L1*Q).
# The operators in this example are purely linear, so we don't worry about Array(µ*L1*Q)[2]

Lₓ = Array(µ*L1*Q)[1] + Array(σ^2/2*L2*Q)[1]
L_bc = I * ρ - Lₓ

# solve the value function
v = L_bc \ π.(x)
