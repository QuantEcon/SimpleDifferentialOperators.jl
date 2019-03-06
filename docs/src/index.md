Installation
==============

To install, run
```julia
] add SimpleDifferentialOperators
```

Note that this requires Julia 1.0 or later.

Usage
==========
Consider solving for `v` from the following equation:
```math
\rho v(x) = f(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

for some constant $\rho, \sigma > 0$ and $\mu \leq 0$. To solve `v` on `M`-size discretized grids, one can run the following code:
```julia
# import LinearAlgebra package (for diagonal and identity matrices)
using LinearAlgebra 
# setup 
f(x) = x^2 
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of grid
x = range(0.0, 1.0, length = M) # grid

# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to reflecting barriers at 0 and 1
bc = (Reflecting(), Reflecting())
L = I * ρ - μ*L₁₋(x, bc) - σ^2 / 2 * L₂(x, bc)
## solve the value function
v_bc = L \ f.(x) 
```

Note that the code above uses differential operators with reflecting boundary conditions applied. 
One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute `v`:
```julia
# operators without boundary conditions, adding extra two rows for boundary conditions
## differential operators on extended nodes
A = μ*L̄₁₋(x) + σ^2 / 2 * L̄₂(x)
## matrix for boundary conditions
B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]]) 
## stack them together
L = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * 0.05 - A); B] 
## solve the value function with reflecting barrier bc (last two elements)
v_bar = L \ [f.(x); 0.0; 0.0] 
## extract the interior (is identical with `v_bc` above)
v_interior = v_bar[2:end-1] 
```

Examples
==========
One can also deploy upwind schemes when drift variable is not constant. Consider solving for `v` from the following equation:
```math
\rho v(x) = f(x) + \mu(x) \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

for the diffusion process
```math
dx = \mu(x) dt + \sigma dW
```

for some constant $\rho, \sigma > 0$ and $\mu(x) = -x$. Note that $\mu(x)$ depends on states. The following code will solve `v` using upwind schemes:
```julia
# setup 
f(x) = x^2 
μ(x) = -x # drift depends on state
σ = 1.0
ρ = 0.05
M = 100 # size of grid
x = range(-1.0, 1.0, length = 100)

bc = (Reflecting(), Reflecting())

# Define first order differential operator using upwind scheme
L₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋(x, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊(x, bc)

# Define linear operator using upwind schemes
L = L₁ - σ^2 / 2 * L₂(x,bc)

# solve the value function
v_bc = (I * ρ - L) \ f.(x) 
```