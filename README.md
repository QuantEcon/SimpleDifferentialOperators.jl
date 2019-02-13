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
μ, σ = (-0.1, 0.1)
z_min, z_max, M = (0.0, 1.0, 100)
grid = range(z_min, z_max, length = M)
L_1_plus, L_1_minus, L_2 = reflecting_diffusionoperators(x) # for "reflecting barrier," or Dirichlet, boundary conditions v'(z_min) = v'(z_max) = 0

# Define a stochastic generator for the process
A = μ*L_1_minus + σ^2 / 2 * L_2 # use L_1_minus because μ < 0  
```

## Documentation

To install, run `] add SimpleDifferentialOperators` on Julia 1.0+.

For more usage information, see the docs badge above.

The [derivation.pdf](docs/tex/derivation.pdf) file in this repository defines and justifies the objects we return.

## Troubleshooting
