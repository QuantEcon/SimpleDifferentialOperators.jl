# SimpleDifferentialOperators

[![Build Status](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl.svg?branch=master)](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl)
[![Codecov](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl)

[![](https://img.shields.io/badge/docs-blue.svg)](https://QuantEcon.github.io/SimpleDifferentialOperators.jl/latest)

## Overview
This is a package to return discretized differential operators of first order and second order under reflecting ("Dirichlet") and mixed (["Robin"](https://en.wikipedia.org/wiki/Robin_boundary_condition)) boundary conditions.

### Example

Consider constructing the corresponding infinitesimal generator for the following stochastic differential equation:
```math
d x_t = μ d_t + σ dW_t  
```
with some constant `μ` and `σ >= 0`, and `W_t` Brownian Motion

```
using SimpleDifferentialOperators
μ = -0.1 # constant negative drift
σ = 0.1
x = range(0.0, 1.0, length = 200) # uniform grid on [0.0, 1.0]

# uses discretized operators for reflecting/Dirichlet boundary conditions, v'(0) = v'(1) = 0
L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())

# discretized generator for the SDE, imposing the boundary conditions
L = μ*L_1_minus + σ^2 / 2 * L_2 # use L_1_minus because μ < 0  
```

## Documentation

To install, run `] add SimpleDifferentialOperators` on Julia 1.1+.

For more usage information, see the docs badge above.

Derivation can be found [here](https://github.com/ubcecon/computing_and_datascience/blob/master/continuous_time_methods/notes/discretized-differential-operator-derivation.tex).

If you want to build the docs locally (say, for contributions), you can just cd to the `docs/` directory and run `julia --project=Project.toml make.jl` (make sure `Documenter.jl` is installed). This will create/populate the `docs/build` directory.

## Troubleshooting

* As a reminder, the package requires **Julia 1.1 or later.**

* If you discover a bug in the code or math, please file an issue in this repo with the label "bug."

* The same holds for feature requests, with the appropriate label.
