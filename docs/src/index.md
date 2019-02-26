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
# setup 
f(x) = x^2 
μ = -0.1 # constant negative drift
σ = 0.1
M = 100 # size of grid
x = range(0.0, 1.0, length = M) # grid

# operators with reflecting boundary conditions
L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())
A = μ*L_1_minus + σ^2 / 2 * L_2 
## solve the value function
v_bc = (I * 0.05 - A) \ f.(x) 
```

Note that the code above uses differential operators with reflecting boundary conditions applied. 
One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute `v`:
```julia
# operators without boundary conditions, adding extra two rows for boundary conditions
## operators on interior nodes
L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x) 
## differential operators on interior nodes
L = μ*L_1_minus + σ^2 / 2 * L_2 
## matrix for boundary conditions
B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]]) 
## stack them together
A = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * 0.05 - L); B] 
## solve the value function with reflecting barrier bc (last two elements)
v_bar = A \ [f.(x); 0.0; 0.0] 
## extract the interior (is identical with `v_bc` above)
v_interior = v_bar[2:end-1] 
```

Examples
==========
One can also deploy upwind schemes when drift variable is not constant. Consider solving for `v` from the following equation:
```math
\rho v(x) = f(x) + \mu(x) \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

for some constant $\rho, \sigma > 0$ and $\mu(x) = -x$. Note that $\mu(x)$ depends on states. The following code will solve `v` using upwind schemes:
```julia
# setup 
f(x) = x^2 
μ(x) = -x # drift depends on state
σ = 0.1
M = 100 # size of grid
x = range(-1.0, 1.0, length = M) # grid
## M-vector of drifts stacked according to the states
μs = μ.(x)

# operators with reflecting boundary conditions
L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())

# Define first order differential operator using upwind scheme
L_1_upwind = (μs .<= 0) .* L_1_minus + (μs .> 0) .* L_1_plus
# Define linear operator using upwind schemes
A = μ.(x) .* L_1_upwind + σ^2 / 2 * L_2 
# solve the value function
v_bc = (I * 0.05 - A) \ f.(x) 
```