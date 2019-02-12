Installation
=====
To install, run
```julia
] dev https://github.com/QuantEcon/SimpleDifferentialOperators.jl
```

Usage 
=====
`robin_diffusionoperators(x, ξ)` takes a vector or range of grids `x` and a coefficient `ξ` to be used for boundary conditions

$$
\begin{aligned}
    \xi v(0) + \nabla_z v(0 ) &= 0 \\
    \xi v(\bar{z}) + \nabla_z v(\bar{z}) &= 0
\end{aligned}
$$
 
and returns a tuple of discretized differential operators in matrix `(L_1_minus, L_1_plus, L_2)` that represent the corresponding forward-difference first-order derivative operator, backward-difference first-order derivative operator, and the centralized second-order derivative operator, respectively. 

For instance, to compute the first-order derivative of scalar-valued function $f$ on given grids $x$ with the boundary condition coefficient $ξ$, one can run the following lines:

```julia
L_1_plus = robin_diffusionoperator(x, ξ) # compute discretized differential operator
f_prime_fd = L_1_plus * f.(x) # first order derivative using forward difference
```

When `ξ` is zero, `reflecting_diffusionoperators(x)` can be used instead.


Derivation
=====
More details on derivation can be found [here](https://github.com/ubcecon/computing_and_datascience/blob/master/continuous_time_methods/notes/differential-operator-on-irregular-grids.tex).

Setup
-------------

-   Define an irregular grid $\{z_i\}_{i=1}^P$ with $z_1 = 0$ and
    $z_P = \bar{z}$ is a "large" number. Denote the grid with the
    variable name, i.e. $z \equiv \{z_i\}_{i=1}^P$.

-   Denote the distance between the grid points as the *backward*
    difference:

    $$
    \begin{aligned}
        \Delta_{i,-} &\equiv z_i - z_{i-1},\, \text{for } i = 2,\ldots, P\\
        \Delta_{i,+} &\equiv z_{i+1} - z_i,\, \text{for } i = 1,\ldots, P-1
    \end{aligned}
    $$

-   Assume $\Delta_{1, -} = \Delta_{1, +}$ and
    $\Delta_{P, +} = \Delta_{P, -}$, due to ghost points, $z_0$ and
    $z_{P+1}$ on both boundaries. (i.e., the distance to the ghost nodes
    are the same as the distance to the closest nodes). Then define the
    vector of backward and forward first differences as 

    $$
        \Delta_{-} \equiv \begin{bmatrix} z_2 - z_1 \\
        \text{diff}(z)
        \end{bmatrix}
    $$

    $$
        \Delta_{+} \equiv \begin{bmatrix} \text{diff}(z)\\
        z_P - z_{P-1}
        \end{bmatrix}
    $$
    

-   Reflecting barrier conditions: 

    $$
    \begin{aligned}
        \xi v(0) + \nabla_z v(0 ) &= 0 & (1) \\
        \xi v(\bar{z}) + \nabla_z v(\bar{z}) &= 0 & (2)
    \end{aligned}
    $$
    
Let $L_1^{-}$ be the discretized backward first differences and $L_2$
be the discretized central differences subject to the Neumann boundary
conditions such that $L_1^{-} v(t)$ and $L_2 v(t)$ represent the first
and second derivatives of $v(z)$ respectively at $t$. For second
derivatives, we use the following numerical scheme from [Achdou et al.
(2017)](http://www.princeton.edu/~moll/HACT.pdf):

$$v''(z_i) \approx \dfrac{ \Delta_{i,-} v( z_i + \Delta_{i,+}) - (\Delta_{i,+} + \Delta_{i,-}) v( z_i ) + \Delta_{i,+} v( z_i - \Delta_{i,-})}{\frac{1}{2}(\Delta_{i,+} + \Delta_{i,-}) \Delta_{i,+} \Delta_{i,-} }, \text{for } i = 1, \ldots, P$$

Regular grids
-------------

Suppose that the grids are regular, i.e., elements of $\text{diff}(z)$
are all identical with $\Delta$ for some $\Delta > 0$.

Using the backward first-order difference, the first boundary condition $(1)$ can be
alternatively represented as 

$$
\begin{aligned}
\dfrac{v(0) - v(-\Delta)}{\Delta} &= - \xi v(0)
\end{aligned}
$$

Similarly, using discretized central differences of second orders,
the first boundary condition $(1)$ yields

$$
\begin{aligned}
\dfrac{v (\Delta) - 2 v(0) + v(-\Delta)}{\Delta^2} &=   \dfrac{v(\Delta) - v(0)}{\Delta^2} - \dfrac{1}{\Delta}\dfrac{v (0) - v(-\Delta) }{\Delta}  \\
&= \dfrac{v(\Delta) - v(0)}{\Delta^2} + \dfrac{1}{\Delta} \xi v(0)  \\ 
&= \dfrac{1}{\Delta^2}  (- 1 + \Delta \xi) v(0)  + \dfrac{1}{\Delta^2}  v(\Delta)  
\end{aligned}
$$

Similarly, for the secondary boundary condition $(2)$, we have 

$$
\begin{aligned}
\dfrac{v (\bar{z} + \Delta) - 2 v(\bar{z} ) + v(\bar{z} -\Delta)}{\Delta^2} &=   \dfrac{v(\bar{z} - \Delta) - v(\bar{z})}{\Delta^2} + \dfrac{1}{\Delta}\dfrac{ v(\bar{z}+\Delta) - v (\bar{z}) }{\Delta}  \\
&= \dfrac{v(\bar{z} - \Delta) - v(\bar{z})}{\Delta^2}  - \dfrac{1}{\Delta} \xi v(\bar{z})  \\ 
&= \dfrac{1}{\Delta^2}  (- 1 - \Delta \xi) v(\bar{z})  + \dfrac{1}{\Delta^2}  v(\bar{z} - \Delta)  
\end{aligned}
$$
 
Thus, the corresponding $L_1^{-}$ and $L_2$ matrices are defined as

$$
\begin{aligned}
L_1^{-} &\equiv \frac{1}{\Delta}\begin{pmatrix}
1 - (1 + \xi \Delta) &0&0&\dots&0&0&0\\
-1&1&0&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&-1&1&0\\
0&0&0&\cdots&0&-1&1
\end{pmatrix}_{P\times P}\\
L_2 &\equiv \frac{1}{\Delta^2}\begin{pmatrix}
-2 + (1 + \xi\Delta) &1&0&\dots&0&0&0\\
1&-2&1&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&1&-2&1\\
0&0&0&\cdots&0&1&-2 + (1- \xi\Delta)
\end{pmatrix}_{P\times P}
\end{aligned}
$$

Irregular grids
---------------

Using the backward first-order difference, the first boundary condition $(1)$ can be
alternatively represented as 

$$
\begin{aligned}
\dfrac{v(0) - v(-\Delta_{1, -})}{\Delta_{1, -}} &= - \xi v(0)
\end{aligned}
$$

Note that we have assumed that $\Delta_{1,-} = \Delta_{1,+}$ and
$\Delta_{P,+} = \Delta_{P,-}$ for the ghost notes. Using discretized
central differences of second orders, the first boundary condition $(1)$ yields

$$
\begin{aligned}
&\dfrac{\Delta_{1,-} v( \Delta_{1,+}) - (\Delta_{1,+} + \Delta_{1,-}) v( 0 ) + \Delta_{1,+}  v( - \Delta_{1,-})}{\frac{1}{2}(\Delta_{1,+} + \Delta_{1,-}) \Delta_{1,+} \Delta_{1,-} } \\
&=
\dfrac{v (\Delta_{1, +}) - 2 v(0) + v(-\Delta_{1, +})}{\Delta_{1, +}^2} \\ &= \dfrac{v(\Delta_{1, +}) - v(0)}{\Delta_{1, +}^2} - \dfrac{1}{\Delta_{1, +}}\dfrac{v (0) - v(-\Delta_{1, +}) }{\Delta_{1, +}}  \\
&= \dfrac{v(\Delta_{1, +}) - v(0)}{\Delta_{1, +}^2} + \dfrac{1}{\Delta_{i,+}} \xi v(0)  \\ 
&= \dfrac{1}{\Delta_{1, +}^2}  (- 1 + \Delta_{1, +} \xi) v(0)  + \dfrac{1}{\Delta_{1, +}^2}  v(\Delta_{1, +})  
\end{aligned}
$$

Similarly, for the secondary boundary condition $(2)$, we have 

$$
\begin{aligned}
&\dfrac{\Delta_{P,-} v( \bar{z} + \Delta_{P,+}) - (\Delta_{P,+} + \Delta_{P,-}) v(\bar{z} ) + \Delta_{P,+}  v( \bar{z} - \Delta_{P,-})}{\frac{1}{2}(\Delta_{P,+} + \Delta_{P,-}) \Delta_{P,+} \Delta_{P,-} } \\
&=\dfrac{v (\bar{z} + \Delta_{P,-}) - 2 v(\bar{z} ) + v(\bar{z} -\Delta_{P,-})}{\Delta_{P,-}^2} \\
&=   \dfrac{v(\bar{z} - \Delta_{P,-}) - v(\bar{z})}{\Delta_{P,-}^2} + \dfrac{1}{\Delta_{P,-}}\dfrac{ v(\bar{z}+\Delta_{P,-}) - v (\bar{z}) }{\Delta_{P,-}}  \\
&= \dfrac{v(\bar{z} - \Delta_{P,-}) - v(\bar{z})}{\Delta_{P,-}^2}  - \dfrac{1}{\Delta_{P,-}} \xi v(\bar{z})  \\ 
&= \dfrac{1}{\Delta_{P,-}^2}  (- 1 - \Delta_{P,-} \xi) v(\bar{z})  + \dfrac{1}{\Delta_{P,-}^2}  v(\bar{z} - \Delta_{P,-})  
\end{aligned}
$$

Thushe corresponding $L_1^{-}$ and $L_2$ matrices are defined as

$$
\begin{aligned}
L_1^{-} &\equiv \begin{pmatrix}
\Delta^{-1}_{1,-} [1 - (1 + \xi \Delta^{-1}_{1,-})] &0&0&\dots&0&0&0\\
-\Delta_{2,-}^{-1}&\Delta_{2,-}^{-1}&0&\dots&0&0&0\\
\vdots&\vdots&\vdots&\ddots&\vdots&\vdots&\vdots\\
0&0&0&\dots&-\Delta_{P-1,-}^{-1}&\Delta_{P-1,-}^{-1}&0\\
0&0&0&\cdots&0&-\Delta_{P,-}^{-1}&\Delta_{P,-}^{-1}
\end{pmatrix}_{P\times P} \\
L_2 &\equiv \begin{pmatrix}
\Delta_{1,+}^{-2}[-2 + (1+\xi \Delta_{1,+})] &\Delta_{1,+}^{-2}&0&\cdots&0&0&0 \\
\vdots&\ddots&\ddots&\ddots&\ddots&\vdots&\vdots\\
0&\cdots&2(\Delta_{i,+}+\Delta_{i,-})^{-1} \Delta_{i,-}^{-1} &-2\Delta_{i,-}^{-1} \Delta_{i,+}^{-1}  & 2 (\Delta_{i,+}+\Delta_{i,-})^{-1} \Delta_{i,+}^{-1}&\cdots&0 \\
\vdots&\vdots&\vdots&\ddots&\ddots&\ddots&\vdots\\
0&0&0&\cdots&0&\Delta_{P,-}^{-2}&\Delta_{P,-}^{-2} [-2 + (1- \xi\Delta_{P,-})]
\end{pmatrix}_{P\times P}
\end{aligned}
$$
