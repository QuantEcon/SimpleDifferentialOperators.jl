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
v_bc = (I * 0.05 - A) \ f.(x) # solve the value function
```

Note that the code above uses differential operators with reflecting boundary conditions applied. 
One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute `v`:
```julia
# operators without boundary conditions, adding extra two rows for boundary conditions
L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x) # operators on interior nodes
L = μ*L_1_minus + σ^2 / 2 * L_2 # differential operators on interior nodes
B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]]) # matrix for boundary conditions
A = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * 0.05 - L); B] # stack them
v_bar = A \ [f.(x); 0.0; 0.0] # solve the value function with reflecting barrier bc
v_interior = v_bar[2:end-1] # extract the interior (is identical with `v_bc`)
```
