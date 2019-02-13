Installation
==============

To install, run
```julia
] add SimpleDifferentialOperators
```

Note that this requires Julia 1.0 or later.

For the complete derivations of the objects we return, see the [PDF](../tex/derivation.pdf).

Usage
==========

This package provides simple differential operators to discretize the [diffusion problem](www.princeton.edu/~moll/HACTproject/option_simple.pdf):

```math
d x_t = \mu(x_t) d_t + \sigma(x_t)d W_t
```

We provide two methods.

[DEPENDENT ON FINAL METHOD STRUCTURE]
