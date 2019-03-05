# SimpleDifferentialOperators

[![Build Status](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl.svg?branch=master)](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl)
[![Codecov](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://QuantEcon.github.io/SimpleDifferentialOperators.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://QuantEcon.github.io/SimpleDifferentialOperators.jl/dev)

## Overview
This is a package to return discretized differential operators subject to various boundary conditions.  It is intended to be a "simple" stopgap as more advanced implementations (e.g. [DiffEqOperators.jl](https://github.com/JuliaDiffEq/DiffEqOperators.jl/) ) mature.

### Example

Consider constructing the corresponding infinitesimal generator for the following stochastic differential equation:
```math
d x_t = μ d_t + σ dW_t
```
with some constant `μ` and `σ >= 0`, and `W_t` Brownian Motion subject, with reflecting barriers at `x=0` and `x=1`.

If the payoff is in state `x` is `f(x)` and `ρ` is the discount rate, then the bellman equation for the expected present discounted value of payoffs fulfills

![Bellman](https://quicklatex.com/cache3/c8/ql_924229b2e52b5323f8133279ce2e2ac8_l3.png)
<!--  \rho v(x) = f(x) + \partial_x v(x) + \frac{\sigma^2}{2}\partial_{xx}v(x) -->
<!-- \partial_x v(0) = 0,\, \partial_x v(1) = 0 -->
subject to ![BC](https://quicklatex.com/cache3/8e/ql_1183a672e909e5a76851d18016a9c68e_l3.png)


Written in operator form, define the differential operator
<!-- \mathcal{L} \equiv \rho - \mu \partial_x - \frac{\sigma^2}{2}\partial_{xx} -->

![Operator](https://quicklatex.com/cache3/6a/ql_1cf7400708d6f645fbf18917daf5296a_l3.png)

and the Bellman equation can then be written as

![Bellman with Operator](https://quicklatex.com/cache3/96/ql_1df64101bb60cb16eb8b0c759b0de496_l3.png)


This package provides components to discretize these sorts of `L` operators.  Solving

```
using SimpleDifferentialOperators, LinearAlgebra
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

## Documentation

To install, run `] add SimpleDifferentialOperators` on Julia 1.1+.

For more usage information, see the docs badge above.

If you want to build the docs locally (say, for contributions), you can just cd to the `docs/` directory and run `julia --project=Project.toml make.jl` (make sure `Documenter.jl` is installed). This will create/populate the `docs/build` directory.

## Troubleshooting

* As a reminder, the package requires **Julia 1.1 or later.**

* If you discover a bug in the code or math, please file an issue in this repo with the label "bug."

* The same holds for feature requests, with the appropriate label.
