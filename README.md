# SimpleDifferentialOperators

[![Build Status](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl.svg?branch=master)](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl)
[![Codecov](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://QuantEcon.github.io/SimpleDifferentialOperators.jl/stable)

## Overview
This is a package to return discretized differential operators subject to various boundary conditions.  It is intended to be a "simple" stopgap as more advanced implementations (e.g. [DiffEqOperators.jl](https://github.com/JuliaDiffEq/DiffEqOperators.jl/) ) mature.  This package is also not intended to provide a "higher-level" interface for constructing the equations.  See [EconPDEs.jl](https://github.com/matthieugomez/EconPDEs.jl/) for a package intended to make translation of the sorts of equations used in economics more direct.

### Example
#### Bellman equation

Consider constructing the corresponding infinitesimal generator for the following stochastic differential equation:
<!-- d x_t = \mu d_t + \sigma d W_t -->
![SDE](https://quicklatex.com/cache3/ea/ql_cff23b548c368d3e69b54d42c1f626ea_l3.png)

with some constant `μ` and `σ >= 0`, and `W_t` Brownian Motion subject, with reflecting barriers at `x=0` and `x=1`.

If the payoff is in state `x` is `f(x)` and `ρ` is the discount rate, then the bellman equation for the expected present discounted value of payoffs fulfills

![Bellman](https://quicklatex.com/cache3/c8/ql_924229b2e52b5323f8133279ce2e2ac8_l3.png)
<!--  \rho v(x) = f(x) + \partial_x v(x) + \frac{\sigma^2}{2}\partial_{xx}v(x) -->
<!-- \partial_x v(0) = 0,\, \partial_x v(1) = 0 -->
subject to ![BC](https://quicklatex.com/cache3/8e/ql_1183a672e909e5a76851d18016a9c68e_l3.png)



Written in operator form, define the differential operators
<!-- \mathcal{L} \equiv \rho - \mu \partial_x - \frac{\sigma^2}{2}\partial_{xx} -->

![Operator](https://quicklatex.com/cache3/bb/ql_b0c2d466a06eb9093cefcfe9bd14dcbb_l3.png)

then the Bellman equation can be written as

![Bellman with Operator](https://quicklatex.com/cache3/96/ql_1df64101bb60cb16eb8b0c759b0de496_l3.png)


This package provides components to discretize differential operators.  To implement directly,

```julia
using SimpleDifferentialOperators, LinearAlgebra
f(x) = x^2 
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
x = range(0.0, 1.0, length = 100)

# discretize L = ρ - μ D_x - σ^2 / 2 D_xx
# subject to reflecting barriers at 0 and 1
bc = (Reflecting(), Reflecting())
L_x = μ*L₁₋(x, bc) - σ^2 / 2 * L₂(x, bc)
L = I * ρ - L_x
## solve the value function
v = L \ f.(x) 
```

#### Kolmogorov forward equation
Likewise, one can also compute the corresponding stationary distribution of `x` above by the Kolmogorov forward equation (KFE), i.e., finding `g(x)` that satisfies

![Stationary distribution from KFE](https://quicklatex.com/cache3/b1/ql_4bcbeb92f8a76079935d5276637c69b1_l3.png)

Written in operator form, define the differential operator

![Operator for KFE](https://quicklatex.com/cache3/4d/ql_610ba8019c30240cf65ed7853af36c4d_l3.png)

and the KFE for the stationary distribution can then be written as

![KFE with Operator](https://quicklatex.com/cache3/74/ql_d6698f4d95a1227d01efafb4afff7474_l3.png)

Note that the operator for the KFE in the original equation is the adjoint operator of the operator for the HJBE, `L_x`, and the correct discretization scheme for `L_x^*` is, analogously, done by taking the transpose of the discretized operator for HJBE, `L_x` (See [Gabaix et al., 2016](https://doi.org/10.3982/ECTA13569) and [Achdou et al., 2017](https://ideas.repec.org/p/nbr/nberwo/23732.html)). Hence, one can find the stationary distribution as follows:

```julia
using Arpack # library for extracting eigenvalues and eigenvectors

# extract eigenvalues and eigenvectors, smallest eigenval in magintute first
λ, ϕ = eigs(transpose(L_x), which = :SM); 
# extract the very first eigenvector (associated with the smallest eigenvalue)
g_ss = real.(ϕ[:,1]);
# normalize it
g_ss = g_ss / sum(g_ss)
```

## Documentation

To install, run `] add SimpleDifferentialOperators` on Julia 1.1+.

For more usage information, see the docs badge above.

If you want to build the docs locally (say, for contributions), you can just cd to the `docs/` directory and run `julia --project=Project.toml make.jl` (make sure `Documenter.jl` is installed). This will create/populate the `docs/build` directory.

## Troubleshooting

* As a reminder, the package requires **Julia 1.1 or later.**

* If you discover a bug in the code or math, please file an issue in this repo with the label "bug."

* The same holds for feature requests, with the appropriate label.
