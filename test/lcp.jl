# VANILLA EXAMPLE
StoppingProblem = (μ_bar = -0.01,
                          σ_bar = 0.01,
                          S_bar = 10.0,
                          γ = 0.5, # u(x) = x^γ
                          ρ = 0.05, # discount rate
                          x_min = 0.01,
                          x_max = 5.0,
                          M = 300) # num of grid points

function LCP_objects(sp)
# setup
  μ_bar, σ_bar, S_bar, γ, ρ, x_min, x_max, M = sp
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

A, B, S, q = LCP_objects(StoppingProblem)
f = z -> B*z + q
n = StoppingProblem.M
lb = zeros(n)
ub = 300*ones(n) # a reasonable guess?
options(convergence_tolerance = 1e-12, output = :no, time_limit = 600) # 10 minute budget
exit_code, sol_z, sol_f = @suppress solveLCP(f, lb, ub);
