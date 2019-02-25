Installation
==============

To install, run
```julia
] add SimpleDifferentialOperators
```

Note that this requires Julia 1.0 or later.

For the complete derivations of the objects we return, see the [PDF](https://github.com/ubcecon/computing_and_datascience/blob/master/continuous_time_methods/notes/discretized-differential-operator-derivation.tex).

Formulas
==========
Under ``M``-length of grids on ``v`` with end points ``z_{\min} < z_{\max}`` with reflecting barrier conditions of

```math
\begin{align}
\underline{\xi} v(z_{\min}) + \nabla v(z_{\min}) &= 0\\
\overline{\xi} v(z_{\max}) + \nabla v(z_{\max}) &= 0
\end{align}
```

We use the following discretization schemes:

```math
L_1^{-} \equiv \frac{1}{\Delta}\begin{pmatrix}
1 - (1 + \underline{\xi} \Delta) &0&0&\dots&0&0&0\\
-1&1&0&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&-1&1&0\\
0&0&0&\cdots&0&-1&1
\end{pmatrix}_{P\times P}
```

```math
L_1^{+} \equiv \frac{1}{\Delta}\begin{pmatrix}
-1&1&0&\dots&0&0&0\\
0&-1&1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&0&-1&1\\
0&0&0&\cdots&0&0&-1+(1-\overline{\xi} \Delta)
\end{pmatrix}_{P\times P}\label{eq:L-1-plus-regular} \\
```

```math
L_2 \equiv \frac{1}{\Delta^2}\begin{pmatrix}
-2 + (1 + \underline{\xi} \Delta) &1&0&\dots&0&0&0\\
1&-2&1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&1&-2&1\\
0&0&0&\cdots&0&1&-2 + (1- \overline{\xi} \Delta)
\end{pmatrix}_{P\times P}\label{eq:L-2-regular}
```

which represent the backward first order, foward first order, and central second order differential operators respectively.

Derivation, including formula for irregular grids, can be found [here](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/blob/master/docs/tex/discretized-differential-operator-derivation.pdf).


Usage
==========

This package provides discretized differential operators of first order and second order under reflecting boundary conditions.

```@autodocs
Modules = [SimpleDifferentialOperators]
Order   = [:function, :type]
```
