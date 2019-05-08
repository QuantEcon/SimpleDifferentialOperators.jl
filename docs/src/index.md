[SimpleDifferentialOperators.jl](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/)
=============

## Installation

To install, run
```julia
] add SimpleDifferentialOperators
```

Note that this requires Julia 1.1 or later.

## Usage

Detailed derivations and more applications can be found [here](generated/discretized-differential-operator-derivation.pdf).

### Solving HJBE with constant drifts
-------------
Consider solving for `v` from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):

```math
\rho v(x) = f(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

for some constant $\rho, \sigma > 0$ and $\mu \leq 0$. To solve `v` under the reflecting barrier conditions $v'(0) = v'(1) = 0$ on `M`-size discretized grids, one can run the following code:

```julia
using LinearAlgebra, SimpleDifferentialOperators
# setup
f(x) = x^2
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
v = L_bc \ f.(x)
```

Note that the interior solution `v` does not the values of $v$ at the boundary, i.e., $v(0)$ and $v(1)$. To extend the interior solution to the boundary points, one can call `extrapolatetoboundary` as follows:

```julia
̄v = extrapolatetoboundary(x̄, v, bc);
```

Here is a complete plot for `v`:

```julia
using Plots
plot(x̄, v̄, lw = 4, label = "v")
```

![plot-hjbe-both-reflecting](assets/plot-hjbe-both-reflecting.png)

Note that the code above uses differential operators on the interior nodes with reflecting boundary conditions applied.
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
v̄ = [L; B] \ [f.(x); b]

# extract the interior (is identical with `v` above)
v =  v̄[2:end-1]
```

### Solving HJBE with absorbing barrier conditions
Instead of having the reflecting barrier conditions on both lower bound and upper bound $v'(0) = v'(1) = 0$ as above, one can impose an absorbing barrier condition as well. To solve `v` under the reflecting barrier conditions $v(0) = S$ (absorbing barrier on lower bound) for some S and $v'(1) = 0$ (reflecting barrier on upper bound), one can construct `B` and `b` for the boundary conditions as follows.

First, consider the case where $S = 0$.

```julia
# define S
S = 0.0

# boundary conditions (i.e. B v̄ = b)
B = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])
b = [S; 0.0];
```

and solve `v`:
```julia
# stack the systems of bellman and boundary conditions, and solve
v̄ = [L; B] \ [f.(x); b]
``` 

Here is a plot for `v`:

```julia
plot(x̄, v̄, lw = 4, label = "v")
```

![plot-hjbe-lb-absorbing-ub-reflecting](assets/plot-hjbe-lb-absorbing-ub-reflecting.png)

Note that this can be alternatively done by constructing the corresponding differential operators on the interior with `Absorbing()` boundary condition:
```julia
# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to reflecting barriers at 0 and 1
bc = (Absorbing(), Reflecting())
Lₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄ , bc)
L_bc = I * ρ - Lₓ

# solve the value function
v = L_bc \ f.(x)
```

In fact, on the interior, they return identical solutions:
```julia
# confirm that they return the identical solution as the one from the stacked system
using Test
@test v ≈ v̄[2:end-1]
```


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
Lₓ = L₁ - σ^2 / 2 * L₂bc(x̄, bc)
L_bc_state_dependent = I * ρ - Lₓ

# solve the value function
v = L_bc_state_dependent \ f.(x)
```

### Finding stationary distribution from the Kolmogorov forward equation (KFE)
-------------
The KFE equation is
$$
\partial_t v(x,t) = -\mu \partial_{x} v(x,t) + \frac{\sigma^2}{2} \partial_{xx} v(x,t)
$$
for $x \in (x_{\min}, x_{\max})$ with the following corresponding reflecting barrier conditions:
```math
\begin{align}
-\mu v(x_{\min}, t) +\frac{\sigma^2}{2} \partial_{x} v(x_{\min}, t) &= 0 \\
-\mu v(x_{\max}, t) +\frac{\sigma^2}{2} \partial_{x} v(x_{\max}, t) &= 0
\end{align}
```

i.e.,

```math
\begin{align}
-\frac{2\mu}{\sigma^2} v(x_{\min}, t) +\partial_{x} v(x_{\min}, t) &= 0 \\
-\frac{2\mu}{\sigma^2} v(x_{\max}, t) +\partial_{x} v(x_{\max}, t) &= 0
\end{align}
```

which gives mixed boundary conditions with $\overline{\xi} = \underline{\xi} = -\frac{2\mu}{\sigma^2}$.

One can compute the stationary distribution of the state `x` above from the corresponding KFE by taking $\partial_{t} g(x,t) = 0$, i.e., solving $g$ from the $L^* g(x) = 0$ where
```math
L^* = - \mu(x) \partial_{x} + \frac{\sigma^2}{2} \partial_{xx}
```

The following code constructs $L^*$:
```julia
# parameter setup
μ = -0.1 # constant negative drift
σ = 0.1
M = 100 # size of grid (interior points)
x_min = 0.0
x_max = 1.0
x̄ = range(x_min, x_max, length = (M+2))

# ξ values for mixed boundary conditions
ξ_lb = ξ_ub = -2μ/σ^2

# define the corresponding mixed boundary conditions
# note that the direction on the lower bound is backward (default is forward)
# as the drift μ is negative.
bc = (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub))

# use SimpleDifferentialOperators.jl to construct the operator on the interior
L_KFE = Array(-μ*L₁₊bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc))
```

One can find the stationary distribution $g$ by solving the following discretized system of equations:

```math
L^* g = 0
```

such that the sum of $g$ is one. This can be found by finding a non-trivial eigenvector `g_ss` for `L_KFE` associated with the eigenvalue of zero:

```julia
using Arpack # library for extracting eigenvalues and eigenvectors

# extract eigenvalues and eigenvectors, smallest eigenval in magintute first
λ, ϕ = eigs(L_KFE, which = :SM); 
# extract the very first eigenvector (associated with the smallest eigenvalue)
g_ss = real.(ϕ[:,1]);
# normalize it
g_ss = g_ss / sum(g_ss)
```

Using `L` from the state-dependent drift example above, this results in the following stationary distribution:

```julia
plot(x, g_ss, lw = 4, label = "g_ss")
```g

![plot-stationary-dist](assets/plot-stationary-dist.png)

Note that the operator for the KFE in the original equation is the adjoint of the operator for infinitesimal generator used in the HJBE, $L$, and the correct discretization scheme for $L^*$ is, analogously, done by taking the transpose of the discretized operator for HJBE, $L$ (See [Gabaix et al., 2016](https://doi.org/10.3982/ECTA13569) and [Achdou et al., 2017](https://ideas.repec.org/p/nbr/nberwo/23732.html)), which has been constructed as `Lₓ` is the HJBE example above. In fact, the discretized $L^*$ and $L^T$ are identical:

```julia
# discretize L = μ D_x + σ^2 / 2 D_xx
# for infinitesimal generators used in the HJBE
# subject to reflecting barrier conditions
bc = (Reflecting(), Reflecting())
Lₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)

@test transpose(Lₓ) == L_KFE
```