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
# SimpleDifferentialOperators setup
function SDO_negative_drift(π, params)
    bc = (Reflecting(), Reflecting())
    Lₓ = params.μ*L₁₋bc(params.x̄, bc) + params.σ^2 / 2 * L₂bc(params.x̄, bc)
    L_bc = I * params.ρ - Lₓ
    v = L_bc \ π.(params.x);
    return (v = v, L_bc = L_bc, Lₓ = Lₓ, L₁₋bc = L₁₋bc(params.x̄, bc), L₂bc = L₂bc(params.x̄, bc))
end
# DiffEqOperators Setup
function DEO_negative_drift(π, params)
    dx = params.x[2] - params.x[1]
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
    v = L_bc \ π.(params.x);
    return (v = v, L_bc = L_bc, Lₓ = Lₓ, L₁₋bc = Array(L1*Q)[1], L₂bc = Array(L2*Q)[1] )
end
@testset "Constant Negative Drifts" begin
    @test DEO_negative_drift(π, p).v ≈ SDO_negative_drift(π, p).v
    @test DEO_negative_drift(π, p).L_bc ≈ SDO_negative_drift(π, p).L_bc
    @test DEO_negative_drift(π, p).Lₓ ≈ SDO_negative_drift(π, p).Lₓ
    @test DEO_negative_drift(π, p).L₁₋bc ≈ SDO_negative_drift(π, p).L₁₋bc
    @test DEO_negative_drift(π, p).L₂bc ≈ SDO_negative_drift(π, p).L₂bc
end
p = params(μ = 0.1);
# SimpleDifferentialOperators setup
function SDO_positive_drift(π, params)
    bc = (Reflecting(), Reflecting())
    Lₓ = params.μ*L₁₊bc(params.x̄, bc) + params.σ^2 / 2 * L₂bc(params.x̄, bc)
    L_bc = I * params.ρ - Lₓ
    v = L_bc \ π.(params.x);
    return (v = v, L_bc = L_bc, Lₓ = Lₓ, L₁₊bc = L₁₊bc(params.x̄, bc), L₂bc = L₂bc(params.x̄, bc))
end
# DiffEqOperators Setup
function DEO_positive_drift(π, params)
    dx = params.x[2] - params.x[1]
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
    v = L_bc \ π.(params.x);
    return (v = v, L_bc = L_bc, Lₓ = Lₓ, L₁₊bc = Array(L1*Q)[1], L₂bc = Array(L2*Q)[1] )
end

@testset "Constant Positive Drifts" begin
    @test DEO_positive_drift(π, p).v ≈ SDO_positive_drift(π, p).v
    @test DEO_positive_drift(π, p).L_bc ≈ SDO_positive_drift(π, p).L_bc
    @test DEO_positive_drift(π, p).Lₓ ≈ SDO_positive_drift(π, p).Lₓ
    @test DEO_positive_drift(π, p).L₁₊bc ≈ SDO_positive_drift(π, p).L₁₊bc
    @test DEO_positive_drift(π, p).L₂bc ≈ SDO_positive_drift(π, p).L₂bc
end

μ(x) = -x
# SimpleDifferentialOperators setup
function SDO_state_dependent_drift(π, μ, params)
    bc = (Reflecting(), Reflecting())
    L₁ = Diagonal(min.(μ.(params.x), 0.0)) * L₁₋bc(params.x̄, bc) + Diagonal(max.(μ.(params.x), 0.0)) * L₁₊bc(params.x̄, bc)
    Lₓ = L₁ - params.σ^2 / 2 * L₂bc(params.x̄, bc)
    L_bc = I * params.ρ - Lₓ
    v = L_bc \ π.(params.x)
    return (v = v, L_bc = L_bc,  Lₓ = Lₓ, L₁ = L₁, L₂bc = L₂bc(params.x̄, bc))
end
# DiffEqOperators Setup
function DEO_state_dependent_drift(π, μ, params)
    dx = params.x[2] - params.x[1]
    # discretize L = ρ - μ D_x - σ^2 / 2 D_xx
    # subject to reflecting barriers at 0 and 1
    L1 = UpwindDifference(1,1,dx,params.M,t->1.0)
    L2 = CenteredDifference(2,2,dx,params.M)
    Q = Neumann0BC(dx, 1)
    # Here Array(A::GhostDerivativeOperator) will return a tuple of the linear part
    # and the affine part of the operator A, hence we index Array(µ*L1*Q).
    # The operators in this example are purely linear, so we don't worry about Array(µ*L1*Q)[2]
    L₁ = Diagonal(min.(μ.(params.x), 0.0)) * Array(L1*Q)[1] + Diagonal(max.(μ.(params.x), 0.0)) * Array(L1*Q)[1]
    Lₓ = L₁ - Array(params.σ^2/2*L2*Q)[1]
    L_bc = I * params.ρ - Lₓ

    # solve the value function
    v = L_bc \ π.(params.x);
    return (v = v, L_bc = L_bc, Lₓ = Lₓ, L₁ = L₁, L₂bc = Array(L2*Q)[1] )
end

@testset "State Dependent Drifts" begin
    @test DEO_state_dependent_drift(π, μ, p).v ≈ SDO_state_dependent_drift(π, μ, p).v
    @test DEO_state_dependent_drift(π, μ, p).L_bc ≈ SDO_state_dependent_drift(π, μ, p).L_bc
    @test DEO_state_dependent_drift(π, μ, p).Lₓ ≈ SDO_state_dependent_drift(π, μ, p).Lₓ
    @test DEO_state_dependent_drift(π, μ, p).L₁ ≈ SDO_state_dependent_drift(π, μ, p).L₁
    @test DEO_state_dependent_drift(π, μ, p).L₂bc ≈ SDO_state_dependent_drift(π, μ, p).L₂bc
end
