# SimpleDifferentialOperators

[![Build Status](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl.svg?branch=master)](https://travis-ci.com/QuantEcon/SimpleDifferentialOperators.jl)
[![Codecov](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/QuantEcon/SimpleDifferentialOperators.jl)

[![](https://img.shields.io/badge/docs-blue.svg)](https://QuantEcon.github.io/SimpleDifferentialOperators.jl/latest)

## Overview

This is a package to return discretized differential operators to solve the (continuous-time) [diffusion problem](www.princeton.edu/~moll/HACTproject/option_simple.pdf):

```math
d x_t = μ(x_t) d_t + σ dW_t  
```

A simple case might be:

```
using SimpleDifferentialOperators
μ = -0.1 # constant negative drift
σ = 0.1
grid = range(0.0, 1.0, length = 200) # uniform grid on [0.0, 1.0]

# get operators for reflecting/Dirichlet boundary conditions, v'(0) = v'(1) = 0
L_1_plus, L_1_minus, L_2 = reflecting_diffusionoperators(grid)

# (discretized) stochastic generator
A = μ*L_1_minus + σ^2 / 2 * L_2 # use L_1_minus because μ < 0  # discretized stochastic generator
```

## Documentation

To install, run `] add SimpleDifferentialOperators` on Julia 1.0+.

For more usage information, see the docs badge above.

The [derivation.pdf](docs/tex/derivation.pdf) file in this repository defines and justifies the objects we return.

## Troubleshooting

* As a reminder, the package requires **Julia 1.0 or later.**

* If you discover a bug in the code or math, please file an issue in this repo with the label "bug."

* The same holds for feature requests, with the appropriate label.
