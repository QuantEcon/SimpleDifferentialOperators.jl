[SimpleDifferentialOperators.jl](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/)
=============

## Installation

To install, run
```julia
] add SimpleDifferentialOperators
```

Note that this requires Julia 1.1 or later.

## Usage

### Solving HJBE with constant drifts
-------------
Consider solving for `v` from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):
```math
\rho v(x) = f(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

for some constant $\rho, \sigma > 0$ and $\mu \leq 0$. To solve `v` under the reflecting barrier conditions $v'(0) = v'(1) = 0$ on `M`-size discretized grids, one can run the following code:

```julia
# import LinearAlgebra package (for diagonal and identity matrices)
using LinearAlgebra
# setup
f(x) = x^2
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of grid (interior points)

x̄ = range(0.0, 1.0, length = (M+2))
x = x̄[2:end-1]

# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to reflecting barriers at 0 and 1
bc = (Reflecting(), Reflecting())
L = I * ρ - μ*L₁₋bc(x̄, bc) - σ^2 / 2 * L₂bc(x̄, bc)
## solve the value function
v = L \ f.(x)
```

Note that the code above uses differential operators with reflecting boundary conditions applied.
One can alternatively use operators on extended nodes (extended operators) and stack them with matrices for boundary conditions to compute `v`:
```julia
# import SparseArrays package (for identity matrix and spzeros)
using SparseArrays

# differential operators on extended nodes
Lₓ = μ*L₁₋(x̄) + σ^2 / 2 * L₂(x̄)

# boundary conditions (i.e. B v̄ = b)
B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])
b = [0.0; 0.0]

# form bellman equation on extension
L = [spzeros(M) ρ*I spzeros(M)] - Lₓ

# stack the systems of bellman and boundary conditions, and solve
v̄ =  [L; B] \ [f.(x); b]

# extract the interior (is identical with `v` above)
v =  v̄[2:end-1]
```


Here is a plot for `v`:

```julia
using Plots
plot(x, v, lw = 4, label = "v")
```

![plot-hjbe-both-reflecting](assets/plot-hjbe-both-reflecting.png)

### Solving HJBE with absorbing barrier conditions
Instead of having the reflecting barrier conditions on both lower bound and upper bound $v'(0) = v'(1) = 0$ as above, one can impose an absorbing barrier condition as well. To solve `v` under the reflecting barrier conditions $v(0) = S$ (absorbing barrier on lower bound) for some S and $v'(1) = 0$ (reflecting barrier on upper bound), one can construct `B` and `b` for the boundary conditions as follows:

```julia
# define S
S = 3.0

# boundary conditions (i.e. B v̄ = b)
B = transpose([[0; 1; zeros(M)] [zeros(M); -1; 1]])
b = [S; 0.0];
```

and solve `v`:
```julia
# stack the systems of bellman and boundary conditions, and solve
v̄ =  [L; B] \ [f.(x); b]

# extract the interior (is identical with `v` above)
v =  v̄[2:end-1]
```

Note that this can be alternatively done by

Here is a plot for `v`:

```julia
plot(x, v, lw = 4, label = "v")
```

![plot-hjbe-lb-absorbing-ub-reflecting](assets/plot-hjbe-lb-absorbing-ub-reflecting.png)

### Solving HJBE with state-dependent drifts
-------------
One can also deploy upwind schemes when drift variable is not constant. Consider solving for `v` from the following Bellman equation:
```math
\rho v(x) = f(x) + \mu(x) \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

associated with the diffusion process
```math
dx = \mu(x) dt + \sigma dW
```

for some constant $\rho, \sigma > 0$ and $\mu(x) = -x$. Note that $\mu(x)$ depends on states. The following code will solve `v` using upwind schemes, with the reflecting barrier conditions $v'(0) = v'(1) = 0$ applied:

```julia
# setup
f(x) = x^2
μ(x) = -x # drift depends on state
σ = 1.0
ρ = 0.05
M = 100 # size of grid

x̄ = range(-1., 1., length = M + 2)
x = interiornodes(x̄) # i.e., x̄[2:end-1]

bc = (Reflecting(), Reflecting())

# Define first order differential operator using upwind scheme
L₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋bc(x̄, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊bc(x̄, bc)

# Define linear operator using upwind schemes
L_x = L₁ - σ^2 / 2 * L₂bc(x̄, bc)
L = I * ρ - L_x

# solve the value function
v = L \ f.(x)
```
