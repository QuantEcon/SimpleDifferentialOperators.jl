
using LinearAlgebra, SimpleDifferentialOperators, DiffEqOperators, Test, Parameters
#parameters
params = @with_kw (
    μ = -0.1, # constant negative drift
    σ = 0.1,
    ρ = 0.05,
    M = 3, # size of grid (interior points)
    x̄ = range(0.0, 1.0, length = (M+2)),
    x = interiornodes(x̄) # i.e., x̄[2:end-1]
)
p = params();
# payoff function
π(x) = x^2
# Below is a basic example from SDO:
function SDO(π, params)
    x = interiornodes(params.x̄) # i.e., x̄[2:end-1]

    # discretize L = ρ - μ D_x - σ^2 / 2 D_xx
    # subject to reflecting barriers at 0 and 1
    bc = (Reflecting(), Reflecting())
    Lₓ = params.μ*L₁₋bc(params.x̄, bc) + params.σ^2 / 2 * L₂bc(params.x̄, bc)
    L_bc = I * params.ρ - Lₓ

    # solve the value function
    v = L_bc \ π.(x);
    return v
end

# The DiffEq approach (as of the current implementation)
function DEO(π, params)
    x = interiornodes(params.x̄)
    dx = x[2] - x[1]
    # discretize L = ρ - μ D_x - σ^2 / 2 D_xx
    # subject to reflecting barriers at 0 and 1
    L1 = UpwindDifference(1,1,dx,params.M,t->1.0)
    L2 = CenteredDifference(2,2,dx,params.M)
    Q = Neumann0BC(dx, 1)
    # Here Array(A::GhostDerivativeOperator) will return a tuple of the linear part
    # and the affine part of the operator A, hence we index Array(µ*L1*Q).
    # The operators in this example are purely linear, so we don't worry about Array(µ*L1*Q)[2]
    Lₓ = Array(params.µ*L1*Q)[1] + Array(params.σ^2/2*L2*Q)[1]
    L_bc = I * params.ρ - Lₓ

    # solve the value function
    v = L_bc \ π.(x);
    return v
end
@testset "Constant Drifts" begin
    @test DEO(π, p) ≈ SDO(π, p)
end
