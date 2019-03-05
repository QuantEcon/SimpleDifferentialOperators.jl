# Setup
StoppingProblem = @with_kw (μ_bar = -0.01,
                          σ_bar = 0.01,
                          S_bar = 10.0,
                          γ = 0.5, # u(x) = x^γ
                          ρ = 0.05, # discount rate
                          x_min = 0.01,
                          x_max = 5.0,
                          M = 15) # num of grid points

function LCP_objects(sp)
# setup
  @unpack μ_bar, σ_bar, S_bar, γ, ρ, x_min, x_max, M = sp
  grid = range(x_min, x_max, length = M)
  μ = μ_bar
  S(x) = S_bar
  u(x) = x^γ
  σ = σ_bar
# construct operator
  A = μ * L₁₋(grid, (Reflecting(), Reflecting())) + σ^2/2 * L₂(grid, (Reflecting(), Reflecting()))
# construct/return LCP objects
  B = ρ*I - A
  S_vec = S.(grid)
  q_vec = -u.(grid) + B*S_vec
  return (A = A, B = B, S = S_vec, q = q_vec)
end

function LCPsolve(sp)
  @unpack A, B, S, q = LCP_objects(sp)
  f = z -> B*z + q
  n = sp.M
  lb = zeros(n)
  ub = 300*ones(n) # a reasonable guess?
  options(convergence_tolerance = 1e-12, output = :no, time_limit = 600) # 10 minute budget
  exit_code, sol_z, sol_f = @suppress solveLCP(f, lb, ub);
end

@testset "Backward Case, Nonzero S" begin
  code, sol, f = LCPsolve(StoppingProblem())
  @test code == :Solved
  @test sol[1] == 0.0
  @test sol[5] ≈ 12.066957816758809
  @test f[1] ≈ 0.39946441597137766
end

@testset "Backward Case, Zero S" begin
  code, sol, f = LCPsolve(StoppingProblem(S_bar = 0.))
  @test code == :Solved
  @test sol[1] ≈ 2.050665004133949
  @test sol[8] ≈ 30.258918534086924
  @test f[1] == 0.0
end
