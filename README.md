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
![SDE](https://render.githubusercontent.com/render/math?math=%24dx_t%20%3D%20%5Cmu%20dt%20%2B%20%5Csigma%20dWt%24)

with some constant `μ` and `σ >= 0`, and `W_t` Brownian Motion subject, with reflecting barriers at `x=0` and `x=1`, i.e., `v'(0) = v'(1) = 0`.

If the payoff is in state `x` is `f(x)` and `ρ` is the discount rate, then the bellman equation for the expected present discounted value of payoffs fulfills

<!--  \rho v(x) = f(x) + \partial_x v(x) + \frac{\sigma^2}{2}\partial_{xx}v(x) -->
![Bellman](https://render.githubusercontent.com/render/math?math=%24%5Crho%20v(x)%20%3D%20f(x)%20%2B%20%5Cpartial_x%20v(x)%20%2B%20%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cpartial_%7Bxx%7Dv(x)%24)
<!-- \partial_x v(0) = 0,\, \partial_x v(1) = 0 -->
subject to ![BC](https://render.githubusercontent.com/render/math?math=%24%5Cpartial_x%20v(0)%20%3D%200%2C%5C%2C%20%5Cpartial_x%20v(1)%20%3D%200%20%24)



Written in operator form, define the differential operators
<!-- \mathcal{L} \equiv \rho - \mu \partial_x - \frac{\sigma^2}{2}\partial_{xx} -->
![Operator](https://render.githubusercontent.com/render/math?math=%5Cmathcal%7BL%7D%20%5Cequiv%20%5Crho%20-%20%5Cmu%20%5Cpartial_x%20-%20%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cpartial_%7Bxx%7D)

then the Bellman equation can be written as
<!-- \mathcal{L} v(x) = f(x) -->
![Bellman with Operator](https://render.githubusercontent.com/render/math?math=%24%5Cmathcal%7BL%7D%20v(x)%20%3D%20f(x)%24)

This package provides components to discretize differential operators.  To implement directly,

```julia
using SimpleDifferentialOperators, LinearAlgebra
f(x) = x^2
μ = -0.1 # constant negative drift
σ = 0.1
ρ = 0.05
M = 100 # size of interior nodes
x̄ = range(0.0, 1.0, length = (M+2)) # extended grid
x = interiornodes(x̄) # interior grid

# discretize L = ρ - μ D_x - σ^2 / 2 D_xx on the interior
# subject to reflecting barriers at 0 and 1
bc = (Reflecting(), Reflecting())
L_x = μ*L₁₋bc(x̄, bc) - σ^2 / 2 * L₂bc(x̄, bc)
L = I * ρ - L_x
## solve the value function on the interior
v = L \ f.(x)
```

To extrapolate the interior solution to the boundary, one can call `extrapolatetoboundary` as follows:

```julia
v̄ = extrapolatetoboundary(x̄, v, bc) 
```

## Documentation

To install, run `] add SimpleDifferentialOperators` on Julia 1.1+.

For more usage information, see the docs badge above.

If you want to build the docs locally (say, for contributions), you can just cd to the `docs/` directory and run `julia --project=Project.toml make.jl` (make sure `Documenter.jl` is installed). This will create/populate the `docs/build` directory.

Detailed derivations and more applications can be found [here](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/releases/download/dev/discretized-differential-operator-derivation.pdf).

## Troubleshooting

* As a reminder, the package requires **Julia 1.1 or later.**

* If you discover a bug in the code or math, please file an issue in this repo with the label "bug."

* The same holds for feature requests, with the appropriate label.
