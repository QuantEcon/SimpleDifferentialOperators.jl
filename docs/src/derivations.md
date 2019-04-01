Derivations
==========
Detailed derivation, including formula for irregular grids, can be found [here](../generated/discretized-differential-operator-derivation.pdf).

Setup
----------
Let $\overline{x} = \{x_i\}_{i=0}^{M+1}$ be an extended grid of discretized $M$-length of grids on $x$ with end points $x_0 = x_{\min}$ and $x_{M+1} = x_{\max}$, and $x = \{x_i\}_{i=1}^M$ be the corresponding interior grid. Also, throughout the section, we consider regular grids, i.e., $x_{i+1} - x_i = \Delta$ for some constant $\Delta > 0$ for all $i = 0, ..., M$. Given a real-valued function $v$, let $v(x)$ be the $M$-length vector whose $i$th element is $v(x_i)$. The goal is to construct a matrix $L$ such that $L v(x)$ represents the first-order or second-order derivative of $v$ on $x$ under some boundary conditions.

Mixed Boundary Values
----------
Under $M$-length of grids on $v$ with end points $x_{\min} < x_{\max}$ with mixed boundary value conditions of

```math
\begin{align}
\underline{\xi} v(x_{\min}) + \nabla v(x_{\min}) &= 0 \\
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
Boundary conditions can be applied manually by using operators on extended grids. Define $v(\overline{x})$ as $(M+2)$-vector whose $i$th element is $\overline {x}_{i-1}$. We can then define the following operators on $\overline{x}$:

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
-1 & 1 & 0 & \dots & 0 & 0 & 0 \\
0 & 0 & 0 & \dots & 0 & -1 & 1\\
\end{bmatrix}_{2 \times (M+2)} \quad 
b = \begin{bmatrix}
0 \\
0
\end{bmatrix}
```


Likewise, for non-homogenous boundary conditions $v'(x_{\min}) = \underline{\xi}$ and $v'(x_{\max}) = \overline{\xi}$, one can use

```math
B = \begin{bmatrix}
-1 & 1 & 0 & \dots & 0 & 0 & 0 \\
0 & 0 & 0 & \dots & 0 & -1 & 1\\
\end{bmatrix}_{2 \times (M+2)} \quad 
b = \begin{bmatrix}
 \underline{\xi} \Delta \\
\overline{\xi} \Delta
\end{bmatrix}
```

Applications
-------------
### Hamilton–Jacobi–Bellman equations (HJBE)
-------------
Consider solving for $v$ from the following optimal control problem
```math
v(x_0) = \max_{ {\{\alpha(t) \} }_{t \geq 0} } \int_{0}^\infty e^{-\rho t} r( x(t), \alpha(t )) dt
```

with the law of motion for the state 
```math
dx = \mu dt + \sigma dW 
```


for some constant $\mu \geq 0$ and $\sigma \geq 0$ with $x(0) = x_0$.

Let $\alpha^*(t)$ be the optimal solution. Suppose that $r$ under $\alpha^*(t)$ can be expressed in terms of state variables, $r^* (x)$. Then, the HJBE yields

```math
\rho v(x) = r^*(x) +  \mu  \partial_{x} v(x) + \dfrac{\sigma^2}{2} \partial_{xx} v(x)
```

In terms of differential operators, one can rewrite the equation as
```math
(\rho - \tilde{L}) v(x) = r^*(x)
```

where 

```math
\tilde{L} = \mu \partial_{x} + (\sigma^2/2) \partial_{xx}
```


By descretizing the space of $x$, one can solve the corresponding system by using discretized operators for $\partial_{x}$ ($L_{1+}$), $\partial_{xx}$ ($L_2$) on some grids of length $M$, $\{x_i\}_{i=1}^M$:

```math
L = \mu L_{1+} + \dfrac{\sigma^2}{2} L_{2}
```

so that $v$ under the optimal plan can be computed by solving the following discretized system of equations:

```math
(\rho I - L) v = r^*
```

where $v$ and $r^*$ are $M$-vectors whose $i$th elements are $v(x_i)$ and $r^*(x_i)$, respectively.

#### Full dynamics of distributions
One can also solve the full PDE in KFE equation, given an initial distribution $g(x, 0)$. After discretization, note that the KFE can be rewritten as

```math
\dot{g}(t) = L^T g(t)
```
where $\dot{g}(t)$ is an $M$-vector whose $i$th element is $\partial_{t} g(x_i, t)$, which can be efficently solved by a number of differential equation solvers available in public, including [DifferentialEquations.jl](http://doi.org/10.5334/jors.151).
