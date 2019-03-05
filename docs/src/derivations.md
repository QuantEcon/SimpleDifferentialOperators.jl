Derivations
==========
Detailed derivation, including formula for irregular grids, can be found [here](../generated/discretized-differential-operator-derivation.pdf).

Mixed Boundary Values
----------
Under $M$-length of grids on $v$ with end points $x_{\min} < x_{\max}$ with mixed boundary value conditions of

```math
\begin{align}
\underline{\xi} v(x_{\min}) + \nabla v(x_{\min}) &= 0\\
\overline{\xi} v(x_{\max}) + \nabla v(x_{\max}) &= 0
\end{align}
```

Note that this can be extended to reflecting boundary conditions by assigning $\underline{\xi} = 0$ and $\overline{\xi} = 0$. We use the following discretization schemes:

```math
L_{1-} \equiv \frac{1}{\Delta}\begin{pmatrix}
1 - (1 + \underline{\xi} \Delta) &0&0&\dots&0&0&0\\
-1&1&0&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&-1&1&0\\
0&0&0&\cdots&0&-1&1
\end{pmatrix}_{M\times M}
```

```math
L_{1+} \equiv \frac{1}{\Delta}\begin{pmatrix}
-1&1&0&\dots&0&0&0\\
0&-1&1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&0&-1&1\\
0&0&0&\cdots&0&0&-1+(1-\overline{\xi} \Delta)
\end{pmatrix}_{M\times M}\label{eq:L-1-plus-regular} \\
```

```math
L_2 \equiv \frac{1}{\Delta^2}\begin{pmatrix}
-2 + (1 + \underline{\xi} \Delta) &1&0&\dots&0&0&0\\
1&-2&1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&1&-2&1\\
0&0&0&\cdots&0&1&-2 + (1- \overline{\xi} \Delta)
\end{pmatrix}_{M\times M}\label{eq:L-2-regular}
```

which represent the backward first order, foward first order, and central second order differential operators respectively.
