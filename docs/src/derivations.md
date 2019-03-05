Derivations
==========
Detailed derivation, including formula for irregular grids, can be found [here](../generated/discretized-differential-operator-derivation.pdf).

Setup
----------
Let $\{x_i\}_{i=1}^M$ be a collection of discretized $M$-length of grids on $x$ with end points $x_1 = x_{\min}$ and $x_M = x_{\max}$. Also, throughout the section, we consider regular grids, i.e., $x_{i+1} - x_i = \Delta$ for some constant $\Delta > 0$ for all $i = 1, ..., M-1$. Also, given a real-valued function $v$, let $v(x)$ be the $M$-length vector whose $i$th element is $v(x_i)$. The goal is to construct a matrix $L$ such that $L v(x)$ represents the first-order or second-order derivative of $v$ on $x$ under some boundary conditions.

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

Applying boundary conditions with operators on extended grids
----------
Boundary conditions can be applied manually by using operators on extended grids. This can be done by first extending $x = \{x_i\}_{i=1}^M$ to $ \overline{x} = \{x_i\}_{i=0}^{M+1}$ where $x_{i+1} - x_i = \Delta$ for some constant $\Delta > 0$ for all $i = 0, ..., M$. We call $x_0$ and $x_{M+1}$, extra nodes just before and after $x_{\min}$ and $x_{\max}$, as ghost nodes. Likewise, define $v(\overline{x})$ as $(M+2)$-vector whose $i$th element is $\overline {x}_i$. We can then define the following operators on $\overline{x}$:

```math
\overline{L}_{1-} \equiv \frac{1}{\Delta}\begin{pmatrix}
-1&1&0&\dots&0&0&0\\
0&-1&1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&1&0&0\\
0&0&0&\cdots&-1&1&0
\end{pmatrix}_{M\times (M+2)}
```

```math
\overline{L}_{1+} \equiv \frac{1}{\Delta}\begin{pmatrix}
0&-1&1&\dots&0&0&0\\
0&0&-1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&-1&1&0\\
0&0&0&\cdots&0&-1&1
\end{pmatrix}_{M\times (M+2)} 
```

```math
\overline{L}_2 \equiv \frac{1}{\Delta^2}\begin{pmatrix}
-1&2&-1&\dots&0&0&0\\
0&-1&2&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&2&-1&0\\
0&0&0&\cdots&-1&2&1
\end{pmatrix}_{M\times (M+2)}
```

Suppose that we want to solve a system $L v({x}) = f(x) $ where $L$ is a linear combination of discretized differential operators for some $f(x)$ that represents the values of a function $f$ on discretized $x$. To solve the system under boundary conditions on $v$, one can construct and solve the following extended system:

```math
\begin{bmatrix}
\overline{L} \\
B
\end{bmatrix} 
v(\overline{x}) = 
\begin{bmatrix}
f(x) \\
b
\end{bmatrix} 
```

with $M_E$ by $(M+2)$ matrix $B$ and $M_E$-length vector $b$ that represent the current boundary conditions, where $M_E$ is the number of boundary conditions to be applied. For instance, to apply reflecting barrier conditions $v'(x_{\min}) = v'(x_{\max}) = 0$, one can use

```math
B = \begin{bmatrix}
1 & -1 & 0 & \dots & 0 & 0 & 0 \\
0 & 0 & 0 & \dots & 0 & -1 & 1\\
\end{bmatrix}_{2 \times (M+2)} \quad 
b = \begin{bmatrix}
0 \\
0
\end{bmatrix}
```


Likewise, for mixed boundary conditions $v'(x_{\min}) = \underline{\xi}$ and $v'(x_{\max}) = \overline{\xi}$, one can use

```math
B = \begin{bmatrix}
1 & -1 & 0 & \dots & 0 & 0 & 0 \\
0 & 0 & 0 & \dots & 0 & -1 & 1\\
\end{bmatrix}_{2 \times (M+2)} \quad 
b = \begin{bmatrix}
 \underline{\xi} \Delta \\
\overline{\xi} \Delta
\end{bmatrix}
```
