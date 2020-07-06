using LinearAlgebra, SimpleDifferentialOperators, DiffEqOperators, Test, Parameters
using Arpack
#parameters
params = @with_kw (
    μ = -0.1, # constant negative drift
    σ = 0.1,
    ρ = 0.05,
    M = 3, # size of grid (interior points)
    x̄ = range(0.0, 1.0, length = (M+2)),
    x = interiornodes(x̄), # i.e., x̄[2:end-1]
    S = 3.0
)
p = params();

#----------------------------------
#Testing For Negative Drifts
#----------------------------------
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
#----------------------------------
#Testing For Positive Drifts
#----------------------------------
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


#----------------------------------
#Testing For State Dependent Drifts
#----------------------------------
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



#----------------------------------
#Testing For Absorbing BC
#----------------------------------


# payoff function
π(x) = x^2
# SimpleDifferentialOperators setup
function SDO_absorbing_bc(π, params)
    bc = (NonhomogeneousAbsorbing(params.S), Reflecting())
    Lₓbc = params.μ*L₁₋bc(params.x̄, bc) + params.σ^2 / 2 * L₂bc(params.x̄ , bc)
    L_bc = I * params.ρ - Lₓbc
    # construct the RHS with affine boundary
    π_star = π.(params.x) + Laffine(L_bc, π.(params.x), bc)
    # solve the value function
    v = L_bc \ π_star
    return (v = v, L_bc = L_bc, Lₓbc = Lₓbc, L₁₋bc = L₁₋bc(params.x̄, bc), L₂bc = L₂bc(params.x̄, bc),
        π_star = π_star, π = π.(params.x))
end


# DiffEqOperators Setup
function DEO_absorbing_bc(π, params)
    dx = params.x[2] - params.x[1]

    L1 = UpwindDifference(1,1,dx,params.M,t->1.0)
    L2 = CenteredDifference(2,2,dx,params.M)
    # RobinBC(l::NTuple{3,T}, r::NTuple{3,T}, dx::T, order = 1)
    # The variables in l are [αl, βl, γl], and correspond to a BC of the form αl*u(0) + βl*u'(0) = γl
    # imposed on the lower index boundary. The variables in r are [αr, βr, γr],
    # and correspond to an analagous boundary on the higher index end.
    l = (1.0, 0.0, p.S)
    r = (0.0, 1.0, 0.0)
    Q = RobinBC(l, r, dx, 1)

    Lₓbc = Array(params.µ*L1*Q)[1] + Array(params.σ^2/2*L2*Q)[1]
    L_bc = I * params.ρ - Lₓbc

    π_star = π.(params.x)
    # I had to do it by handcode, I couldn't find something like Laffine in SDO.
    π_star[1] -= params.S*(params.μ/dx - (params.σ^2/2 )/(dx^2) )

    # solve the value function
    v = L_bc \ π.(params.x);
    return (v = v, L_bc = L_bc, Lₓbc = Lₓbc, L₁₋bc = Array(L1*Q)[1], L₂bc = Array(L2*Q)[1],
        π_star = π_star, π = π.(params.x))
end

p = params(x̄ = range(0.0, 1.0, length = (p.M+2)))

@testset "Absorbing BC" begin
    @test DEO_absorbing_bc(π, p).v ≈ SDO_absorbing_bc(π, p).v
    @test DEO_absorbing_bc(π, p).L_bc ≈ SDO_absorbing_bc(π, p).L_bc
    @test DEO_absorbing_bc(π, p).Lₓbc ≈ SDO_absorbing_bc(π, p).Lₓbc
    @test DEO_absorbing_bc(π, p).L₁₋bc ≈ SDO_absorbing_bc(π, p).L₁₋bc
    @test DEO_absorbing_bc(π, p).L₂bc ≈ SDO_absorbing_bc(π, p).L₂bc
    @test DEO_absorbing_bc(π, p).π_star ≈ SDO_absorbing_bc(π, p).π_star
    @test DEO_absorbing_bc(π, p).π ≈ SDO_absorbing_bc(π, p).π
end



#----------------------------------
#Testing For Solving KFE
#----------------------------------


function SDO_Solve_KFE(params)
    # ξ values for mixed boundary conditions
    ξ_lb = ξ_ub = -2params.μ/params.σ^2
    # define the corresponding mixed boundary conditions
    # note that the direction on the lower bound is backward (default is forward)
    # as the drift μ is negative.
    bc = (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub))

    # use SimpleDifferentialOperators.jl to construct the operator on the interior
    L_KFE_with_drift = Array(-params.μ*L₁₊bc(params.x̄, bc) + params.σ^2 / 2 * L₂bc(params.x̄, bc))
    L_KFE_without = params.σ^2 / 2 * L₂bc(params.x̄, bc)
    return (L_KFE_with_drift = L_KFE_with_drift, L_KFE_without = L_KFE_without)
end


function DEO_Solve_KFE(params)
    dx = params.x[2] - params.x[1]

    L1 = UpwindDifference(1,1,dx,params.M,t->1.0)
    L2 = CenteredDifference(2,2,dx,params.M)

    ξ_lb = ξ_ub = -2params.μ/params.σ^2

    l = (1.0, ξ_lb, 0.0)
    r = (1.0, ξ_ub, 0.0)
    Q = RobinBC(l, r, dx, 1)

    # use SimpleDifferentialOperators.jl to construct the operator on the interior
    L_KFE_with_drift = Array(-params.µ*L1*Q)[1] + Array(params.σ^2/2 *L2*Q)[1]
    L_KFE_without = Array(params.σ^2/2*L2*Q)[1]
    return (L_KFE_with_drift = L_KFE_with_drift, L_KFE_without = L_KFE_without)
end


p = params(x̄ = range(0.0, 1.0, length = (p.M+2)))
@testset "Solving KFE" begin
    @test DEO_Solve_KFE(p).L_KFE_with_drift ≈ SDO_Solve_KFE(p).L_KFE_with_drift
    @test DEO_Solve_KFE(p).L_KFE_without ≈ SDO_Solve_KFE(p).L_KFE_without
end
