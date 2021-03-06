
## Solving HJBE with constant drifts
-------------
Consider solving for `v` from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):

```math
\rho v(x) = \pi(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

for some constant $\rho, \sigma > 0$ and $\mu \leq 0$. To solve `v` under the reflecting barrier conditions $v'(0) = v'(1) = 0$ on `M`-size discretized grids, one can run the following code:

```julia
using LinearAlgebra, SimpleDifferentialOperators
# setup
π(x) = x^2
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
v = L_bc \ π.(x)
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
v̄ =  [L; B] \ [π.(x); b]

# extract the interior (is identical with `v` above)
v =  v̄[2:end-1]
```

## Solving HJBE with absorbing barrier conditions
Instead of having the reflecting barrier conditions on both lower bound and upper bound $v'(0) = v'(1) = 0$ as above, one can impose an absorbing barrier condition as well. Consider solving for `v` from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):

```math
\rho v(x) = \pi(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

but this time under the boundary conditions $v(0) = S$ (absorbing barrier on lower bound) for some S and $v'(1) = 0$ (reflecting barrier on upper bound). The corresponding boundary condition matrices `B` and `b` can be constructed as follows; first, consider the case where $S \neq 0$, which gives a nonhomogenous boundary condition:

```julia
using LinearAlgebra, SimpleDifferentialOperators, SparseArrays
# setup
π(x) = x^2
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of grid (interior points)

x̄ = range(0.0, 1.0, length = (M+2))
x = interiornodes(x̄) # i.e., x̄[2:end-1]

# differential operators on extended nodes
Lₓ = μ*L₁₋(x̄) + σ^2 / 2 * L₂(x̄)

# boundary conditions (i.e. B v̄ = b)
B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])
b = [0.0; 0.0]

# form bellman equation on extension
L = [spzeros(M) ρ*I spzeros(M)] - Lₓ

# define S
S = 3.0

# boundary conditions (i.e. B v̄ = b)
B = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])
b = [S; 0.0];
```

Now solve `v`:
```julia
# stack the systems of bellman and boundary conditions, and solve
v̄ =  [L; B] \ [π.(x); b]
```

Or we can proceed as follows: we can apply one Gaussian elimination step to remove a non-zero element of the first column in $L$, which is $\mu \Delta^{-1} - (\sigma^2/2) \Delta^{-2}$. This can be done by substracting the first row of the stacked system $[L; B]$ by the first row of the system $B = b$ by $\mu \Delta^{-1} - (\sigma^2/2) \Delta^{-2}$. This returns the following identical system:

```math
\begin{bmatrix}
L^* \\
B 
\end{bmatrix}
\begin{bmatrix}
v(x_0) \\
v(x_1) \\ 
\vdots \\
v(x_{M}) \\
v(x_{M+1})
\end{bmatrix}
=
\begin{bmatrix}
π^* \\
b
\end{bmatrix}
```

where

```math
π^* =
\begin{bmatrix}
π(x_1) - S(\mu \Delta^{-1} - (\sigma^2/2) \Delta^{-2})
\\ 
\vdots
\\
π(x_{M})
\end{bmatrix}
```

and $L^*$ is identical with $L$ except the first column, which is a zero column in $L^*$. Note that the second boundary condition on the upper bound is homogeneous, so as in the HJBE with constant drifts example, we can apply the boundary condition by Gaussian elimination and achieve the following system for the solution on interior nodes:

```math
\begin{bmatrix}
L^{B*}
\end{bmatrix}
\begin{bmatrix}
v(x_1)
\\ 
\vdots
\\
v(x_{M}) 
\end{bmatrix}
=
π^* 
```

where $L^{B*}$ is the operator on interior nodes with the boundary conditions applied. In the code, one can solve the system as follows:

```julia
# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to boundary conditions at 0 and 1
bc = (NonhomogeneousAbsorbing(S), Reflecting())
Lₓbc = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄ , bc)
Lbc = I * ρ - Lₓbc

# construct the RHS with affine boundary
π_star = π.(x) + Laffine(Lbc, π.(x), bc)

# solve the value function
v = Lbc \ π_star
```

Here is a plot for `v`:

```julia
using Plots
plot(x̄, v̄, lw = 4, label = "v")
```

![plot-hjbe-lb-absorbing-ub-reflecting](assets/plot-hjbe-lb-absorbing-ub-reflecting.png)

Note that this can be alternatively done by constructing the corresponding differential operators on the interior with `Absorbing()` boundary condition when $S = 0$:
```julia
# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to boundary conditions at 0 and 1
bc = (Absorbing(), Reflecting())
Lₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄ , bc)
L_bc = I * ρ - Lₓ

# solve the value function
v = L_bc \ π.(x)
```

In fact, on the interior, they return identical solutions:
```julia
# define S
S = 0.0

# boundary conditions (i.e. B v̄ = b)
B = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])
b = [S; 0.0];

# stack the systems of bellman and boundary conditions, and solve
v̄ = [L; B] \ [π.(x); b]

# confirm that v returns the identical solution as the one from the stacked system
using Test
@test v ≈ v̄[2:end-1]
```

## Solving HJBE with absorbing barrier conditions on interior nodes
Like the example above, consider the following HJBE with boundary conditions $v(0) = 0$ (absorbing barrier on lower bound) for some S and $v'(1) = 0$ (reflecting barrier on upper bound):

```math
\rho v(x) = \pi(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

but this time under the boundary conditions $v(x) = 0$ for all $x \leq x_{M^*}$ for some $0 < M^* < M$ so that absorbing barrier conditions are applied on interior nodes as well. First, one can construct the extension operator `L` with boundary condition matrices `B` and `b` with the boundary conditions $v(0) = 0$ and $v'(1) = 0$ as follows:

```julia
using LinearAlgebra, SimpleDifferentialOperators, SparseArrays
# setup
π(x) = x^2
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of grid (interior points)

x̄ = range(0.0, 1.0, length = (M+2))
x = interiornodes(x̄) # i.e., x̄[2:end-1]

# differential operators on extended nodes
Lₓ = μ*L₁₋(x̄) + σ^2 / 2 * L₂(x̄)

# boundary conditions (i.e. B v̄ = b)
B = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])
b = [0.0; 0.0]

# form bellman equation on extension
L = [spzeros(M) ρ*I spzeros(M)] - Lₓ

# define S
S = 0.0

# boundary conditions (i.e. B v̄ = b)
B = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])
b = [S; 0.0];
```

Note that under $v(x) = 0$ for all $x \leq x_{M^*}$, the following equation must be satisfied for all $0 \geq m \leq M^*$:

```math
v(x_m) = 0 & \forall 0 \geq m \leq M^*
```

In the matrix form with discretized extended solution $\overline{x}$, this is equivalent to

```math
\begin{bmatrix}
0 & 1 & 0 & \dots & 0 & 0
0 & 0 & 1 & \ddots & 0 & 0
0 & 0 & 0 & \ddots & 0 & 0
0 & 0 & 0 & \ddots & 1 & 0
0 & 0 & 0 & \dots & 0 & 1
\end{bmatrix}
\begin{bmatrix}
\overline{x}_0 \\
\overline{x}_1 \\
\overline{x}_2 \\
\vdots \\
\overline{x}_{M^*-1}
\overline{x}_{M^*}
\end{bmatrix}
=
\begin{bmatrix}
0 \\ 0 \\ \vdots \\ 0 \\ 0
\end{bmatrix}
```
so that `L` should be modified such that its first $m$th rows are all zero rows except $(m,m+1)$th element, which has a value of one for all $1 \leq m \leq M^*$. Likewise, the expression on the right hand side, `π.(x)` has to be modified so that its first $M^*$ elements are all zero.

```julia
M_loc = 40 # cutoff for absorbing barrier conditions on interior nodes
L_modified = copy(L)
RHS = copy(π.(x))

# modify L and RHS
L_modified[:,1:(1+M_loc)] = zeros(M,1+M_loc)
L_modified[1:M_loc,2:(1+M_loc)] = Diagonal(ones(M_loc))
RHS[1:M_loc] = zeros(M_loc)

# stack the systems of bellman and boundary conditions, and solve
v̄ =  [L_modified; B] \ [RHS; b]
```

And here's a plot:

```julia
using Plots
plot(x̄, v̄, lw = 4, label = "v")
```

![plot-hjbe-absorbing-interior](assets/plot-hjbe-absorbing-interior.png)

Alternatively, one can use `Lbc` with `Laffine` family as follows:

```julia
# define boundary conditions
bc = (Absorbing(M_loc), Reflecting())

# differential operators on interior nodes, boundary conditions applied
Lₓbc = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)

# form bellman equation on extension
Lbc = ρ*I - Lₓbc

# construct the RHS with affine boundary
π_star = π.(x) + Laffine(Lbc, π.(x), bc)

# solve the value function
v = Lbc \ π_star

# confirm that v returns the identical solution as the one from the stacked system
using Test
@test v ≈ v̄[2:end-1]
```

## Solving HJBE with jump diffusion
Consider the jump process added to the HJBE with some intensity $\lambda \geq 0$: 

```math
\rho v(x) = \pi(x) + \mu \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x) + \lambda \left[ v(x + J(x) ) - v(x) \right]
```

where the jump process is defined by the jump magnitude defined by $J(x_i)$. In `SimpleDifferentialOperators.jl`, the jump process can be defined as follows:

```julia
# uniform jump
jumpf(x_i) = -0.01
jumpprocess = JumpProcess(x̄, jumpf)
```

Note that, the corresponding indices for destinations will be determined by the nearest neighbor as the domain has to be discretized accordingly. Alternatively, if the jump magnitude is uniform across all cohorts, one can forward the uniform jump magnitude as follows:

```julia
# use the fact that the jump magnitude is uniform across all nodes
jumpprocess = JumpProcess(x̄, -0.01)
```

One can define a jump process manually by providing jump magnitudes in indices as well. If a jump process is defined by the indices on a discretized domain, incurring jumps from $v(x_i)$ to $v(x_{i-1})$ for all $i$ in $2 \leq i \leq M$, one can construct a jump process as follows:

```julia
# length of nodes on the interior
M = length(interiornodes(x̄))
# vector of jumps; ith element represents the jump from ith node in the interior
jumps = -ones(M)
# define jump process 
jumpprocess =  JumpProcess(x̄, jumps)
```

Alternatively, one can define an identical jump process with ease if the jump maginitude in index is uniform across all nodes:

```julia
# use the fact that the jump size is uniform across all nodes
jumpprocess = JumpProcess(x̄, -1)
```

Then one can define the corresponding discretized operator $L_n$ and solve value functions as follows:

```julia
# define jump intensity
λ = 0.6

# construct discretized operators on interior nodes with the bc
Lₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄ , bc) + λ * Lₙbc(x̄, bc, jumpprocess) 
L_bc = I * ρ - Lₓ

# solve the value function
v = L_bc \ π.(x)
```

## Solving HJBE with jump diffusions and Markov chains
Suppose we are asked to solve HJBE with two states ($N=2$) where for each $i$th state with the corresponding differential operator $L_i$ under different payoff functions $\pi_i$ and drifts $\mu_i$, there is a transition intensity of $q_{ij}$ to have state $j$ assigned.

```julia
# setup
# payoff functions
π_1(x) = x^2
π_2(x) = (x-0.01)^2

# constant negative drifts
μ_1 = -0.1
μ_2 = -0.15
λ = 0.6
σ = 0.1
ρ = 0.05
M = 100 # size of grid (interior points)

x̄ = range(0.0, 1.0, length = (M+2))
x = interiornodes(x̄) # i.e., x̄[2:end-1]
```

Let the HJBE in the first state have a jump process $J$ associated while the one for the second state does not. Then we have the following system of differential equations; note that we have $q_{ii} = -q_{ij}$ for $i \neq j$:

```math
\begin{align}
\rho v_1 (x) &= \pi_1(x) + \mu_1 \partial_x v_1(x) + \frac{\sigma^2}{2} \partial_{xx} v_1(x) + \lambda \left[ v_1(x + J(x) ) - v_1(x) \right] + q_{12} [ v_2(x) - v_1(x) ] \\
\rho v_2 (x) &= \pi_2(x) + \mu_2 \partial_x v_2(x) + \frac{\sigma^2}{2} \partial_{xx} v_2(x) + q_{21} [ v_1(x) - v_2(x) ]
\end{align}
```

First, construct $L_1$ and $L_2$ ignoring the Markov chain for transition between the two states; assume that both states have reflecting boundary conditions applied.

```julia
# construct the jump process for the operator in state 1
jumpprocess1 = JumpProcess(x̄, -0.01)

# construct the differential operators for both states
# subject to reflecting barriers at 0 and 1
bc = (Reflecting(), Reflecting())
L_1ₓ = μ_1*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc) + λ * Lₙbc(x̄, bc, jumpprocess1)
L_2ₓ = μ_2*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)
L_1_bc = I * ρ - L_1ₓ
L_2_bc = I * ρ - L_2ₓ
```

Then construct an intensity matrix $Q$, whose $(i,j)$th element represents $q_{ij}$:

```julia
# define intensity matrix for transition
Q = [-0.01 0.01; 0.02 -0.02]
```

Using the discretized operators `L_1_bc` and `L_2_bc` on interior nodes `x` with boundary conditions `bc` applied, one can construct the joint operator `L_bc` with the intensity matrix `Q` as follows:

```julia
# define the corresponding joint operator
L_bc = jointoperator_bc((L_1_bc, L_2_bc), Q)
```

Construct a vector of payoff functions $\pi_1$ and $\pi_2$ stacked together and solve the system using the joint operator constructed above:
```julia
# solve the system
v = L_bc \ [π_1.(x); π_2.(x)]
```

Note that the first $M$ elements represent the discretized solution for $v_1$ and the last $M$ elements represent the one for $v_2$:

```julia
# extract the solution for each state
v_1 = v[1:M]
v_2 = v[(M+1):end]

# plot v_1 and v_2
plot(x, [v_1, v_2], lw = 4, label = ["v_1", "v_2"])
```

![plot-hjbe-two-states](assets/plot-hjbe-two-states.png)

## Solving HJBE with state-dependent drifts
-------------
One can also deploy upwind schemes when drift variable is not constant. Consider solving for `v` from the following Bellman equation:
```math
\rho v(x) = π(x) + \mu(x) \partial_x v(x) + \frac{\sigma^2}{2} \partial_{xx} v(x)
```

associated with the diffusion process
```math
dx = \mu(x) dt + \sigma dW
```

for some constant $\rho, \sigma > 0$ and $\mu(x) = -x$. Note that $\mu(x)$ depends on states. The following code will solve `v` using upwind schemes, with the reflecting barrier conditions $v'(0) = v'(1) = 0$ applied:

```julia
# setup
π(x) = x^2
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
v = L_bc_state_dependent \ π.(x)
```

## Finding stationary distribution from the Kolmogorov forward equation (KFE)
-------------
The KFE equation is
```math
\partial_t f(x,t) = -\mu \partial_{x} f(x,t) + \frac{\sigma^2}{2} \partial_{xx} f(x,t)
```

for $x \in (x_{\min}, x_{\max})$ with the following corresponding reflecting barrier conditions:
```math
\begin{align}
-\mu f(x_{\min}, t) +\frac{\sigma^2}{2} \partial_{x} f(x_{\min}, t) &= 0 \\
-\mu f(x_{\max}, t) +\frac{\sigma^2}{2} \partial_{x} f(x_{\max}, t) &= 0
\end{align}
```

i.e.,

```math
\begin{align}
-\frac{2\mu}{\sigma^2} f(x_{\min}, t) +\partial_{x} f(x_{\min}, t) &= 0 \\
-\frac{2\mu}{\sigma^2} f(x_{\max}, t) +\partial_{x} f(x_{\max}, t) &= 0
\end{align}
```

which gives mixed boundary conditions with $\overline{\xi} = \underline{\xi} = -\frac{2\mu}{\sigma^2}$.

One can compute the stationary distribution of the state `x` above from the corresponding KFE by taking $\partial_{t} f(x,t) = 0$, i.e., solving $f$ from the $L^* f(x) = 0$ where
```math
L^* = - \mu \partial_{x} + \frac{\sigma^2}{2} \partial_{xx}
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

One can find the stationary distribution $f$ by solving the following discretized system of equations:

```math
L^* f = 0
```

such that the sum of $f$ is one. This can be found by finding a non-trivial eigenvector `f_ss` for `L_KFE` associated with the eigenvalue of zero:

```julia
using Arpack # library for extracting eigenvalues and eigenvectors

# extract eigenvalues and eigenvectors, smallest eigenval in magintute first
λ, ϕ = eigs(L_KFE, which = :SM); 
# extract the very first eigenvector (associated with the smallest eigenvalue)
f_ss = real.(ϕ[:,1]);
# normalize it
f_ss = f_ss / sum(f_ss)
```

Using `L` from the state-dependent drift example above, this results in the following stationary distribution:

```julia
plot(x, f_ss, lw = 4, label = "f_ss")
```

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
