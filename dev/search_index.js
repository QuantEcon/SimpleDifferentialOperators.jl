var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "﻿Installation",
    "title": "﻿Installation",
    "category": "page",
    "text": ""
},

{
    "location": "#Installation-1",
    "page": "﻿Installation",
    "title": "﻿Installation",
    "category": "section",
    "text": "To install, run] add SimpleDifferentialOperatorsNote that this requires Julia 1.0 or later."
},

{
    "location": "#Usage-1",
    "page": "﻿Installation",
    "title": "Usage",
    "category": "section",
    "text": "Consider solving for v from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):rho v(x) = f(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v on M-size discretized grids, one can run the following code:# import LinearAlgebra package (for diagonal and identity matrices)\nusing LinearAlgebra \n# setup \nf(x) = x^2 \nμ = -0.1 # constant negative drift\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid\nx = range(0.0, 1.0, length = M) # grid\n\n# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nL = I * ρ - μ*L₁₋(x, bc) - σ^2 / 2 * L₂(x, bc)\n## solve the value function\nv_bc = L \\ f.(x) Note that the code above uses differential operators with reflecting boundary conditions applied.  One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute v:# operators without boundary conditions, adding extra two rows for boundary conditions\n## differential operators on extended nodes\nA = μ*L̄₁₋(x) + σ^2 / 2 * L̄₂(x)\n## matrix for boundary conditions\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]]) \n## stack them together\nL = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * 0.05 - A); B] \n## solve the value function with reflecting barrier bc (last two elements)\nv_bar = L \\ [f.(x); 0.0; 0.0] \n## extract the interior (is identical with `v_bc` above)\nv_interior = v_bar[2:end-1] "
},

{
    "location": "#Examples-1",
    "page": "﻿Installation",
    "title": "Examples",
    "category": "section",
    "text": ""
},

{
    "location": "#Solving-HJBE-with-state-dependent-drift-variables-1",
    "page": "﻿Installation",
    "title": "Solving HJBE with state-dependent drift variables",
    "category": "section",
    "text": "One can also deploy upwind schemes when drift variable is not constant. Consider solving for v from the following Bellman equation:rho v(x) = f(x) + mu(x) partial_x v(x) + fracsigma^22 partial_xx v(x)associated with the diffusion processdx = mu(x) dt + sigma dWfor some constant rho sigma  0 and mu(x) = -x. Note that mu(x) depends on states. The following code will solve v using upwind schemes:# setup \nf(x) = x^2 \nμ(x) = -x # drift depends on state\nσ = 1.0\nρ = 0.05\nM = 100 # size of grid\nx = range(-1.0, 1.0, length = 100)\n\nbc = (Reflecting(), Reflecting())\n\n# Define first order differential operator using upwind scheme\nL₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋(x, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊(x, bc)\n\n# Define linear operator using upwind schemes\nL = L₁ - σ^2 / 2 * L₂(x,bc)\n\n# solve the value function\nv_bc = (I * ρ - L) \\ f.(x) "
},

{
    "location": "#Finding-stationary-distribution-from-the-Kolmogorov-forward-equation-(KFE)-1",
    "page": "﻿Installation",
    "title": "Finding stationary distribution from the Kolmogorov forward equation (KFE)",
    "category": "section",
    "text": "One can also compute the stationary distribution of the state x above from the corresponding KFE:partial_t g(xt) = - mu(x) partial_x g(x t) + fracsigma^22 partial_xx g(xt)by taking partial_t g(xt) = 0, i.e., solving g from the L^* g(x) = 0 whereL^* = - mu(x) partial_x + fracsigma^22 partial_xxBy descretizing the space of x, one can solve the corresponding system by using discretized operators for L^*. Note that the operator for the KFE in the original equation is the adjoint operator of the operator for the HJBE, L, and the correct discretization scheme for L^* is, analogously, done by taking the transpose of the discretized operator for HJBE, L (See Gabaix et al., 2016). Hence, one can find the stationary distribution by solving the following discretized system of equations:L^T g = 0such that the sum of g is one. This can be found by finding a non-trivial eigenvector for L^T  associated with the eigenvalue of zero:using Arpack # library for extracting eigenvalues and eigenvectors\n\n# extract eigenvalues and eigenvectors, smallest eigenval in magintute first\nλ, ϕ = eigs(transpose(L), which = :SM); \n# extract the very first eigenvector (associated with the smallest eigenvalue)\ng_ss = real.(ϕ[:,1]);\n# normalize it\ng_ss = g_ss / sum(g_ss)"
},

{
    "location": "notebooks/#",
    "page": "Notebooks",
    "title": "Notebooks",
    "category": "page",
    "text": ""
},

{
    "location": "notebooks/#Notebooks-1",
    "page": "Notebooks",
    "title": "Notebooks",
    "category": "section",
    "text": "We prepared a gallery of notebooks to showcase the applicability of the package to different economic problems.Optimal Stopping Problems. You can also download the notebook."
},

{
    "location": "derivations/#",
    "page": "Derivations",
    "title": "Derivations",
    "category": "page",
    "text": ""
},

{
    "location": "derivations/#Derivations-1",
    "page": "Derivations",
    "title": "Derivations",
    "category": "section",
    "text": "Detailed derivation, including formula for irregular grids, can be found here."
},

{
    "location": "derivations/#Setup-1",
    "page": "Derivations",
    "title": "Setup",
    "category": "section",
    "text": "Let x_i_i=1^M be a collection of discretized M-length of grids on x with end points x_1 = x_min and x_M = x_max. Also, throughout the section, we consider regular grids, i.e., x_i+1 - x_i = Delta for some constant Delta  0 for all i = 1  M-1. Also, given a real-valued function v, let v(x) be the M-length vector whose ith element is v(x_i). The goal is to construct a matrix L such that L v(x) represents the first-order or second-order derivative of v on x under some boundary conditions."
},

{
    "location": "derivations/#Mixed-Boundary-Values-1",
    "page": "Derivations",
    "title": "Mixed Boundary Values",
    "category": "section",
    "text": "Under M-length of grids on v with end points x_min  x_max with mixed boundary value conditions ofmath \\underline{\\xi} v(x{\\min}) + \\nabla v(x{\\min}) &= 0\\\n\\overline{\\xi} v(x{\\max}) + \\nabla v(x{\\max}) &= 0Note that this can be extended to reflecting boundary conditions by assigning underlinexi = 0 and overlinexi = 0. We use the following discretization schemes:L_1- equiv frac1Deltabeginpmatrix\n1 - (1 + underlinexi Delta) 00dots000\n-110dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Mtimes ML_1+ equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots0-11\n000cdots00-1+(1-overlinexi Delta)\nendpmatrix_Mtimes MlabeleqL-1-plus-regular L_2 equiv frac1Delta^2beginpmatrix\n-2 + (1 + underlinexi Delta) 10dots000\n1-21dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots1-21\n000cdots01-2 + (1- overlinexi Delta)\nendpmatrix_Mtimes MlabeleqL-2-regularwhich represent the backward first order, foward first order, and central second order differential operators respectively."
},

{
    "location": "derivations/#Applying-boundary-conditions-with-operators-on-extended-grids-1",
    "page": "Derivations",
    "title": "Applying boundary conditions with operators on extended grids",
    "category": "section",
    "text": "Boundary conditions can be applied manually by using operators on extended grids. This can be done by first extending x = x_i_i=1^M to $ \\overline{x} = {xi}{i=0}^{M+1}$ where x_i+1 - x_i = Delta for some constant Delta  0 for all i = 0  M. We call x_0 and x_M+1, extra nodes just before and after x_min and x_max, as ghost nodes. Likewise, define v(overlinex) as (M+2)-vector whose ith element is overline x_i. We can then define the following operators on overlinex:overlineL_1- equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots100\n000cdots-110\nendpmatrix_Mtimes (M+2)overlineL_1+ equiv frac1Deltabeginpmatrix\n0-11dots000\n00-1dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Mtimes (M+2) overlineL_2 equiv frac1Delta^2beginpmatrix\n-12-1dots000\n0-12dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots2-10\n000cdots-121\nendpmatrix_Mtimes (M+2)Suppose that we want to solve a system L v(x) = f(x)  where L is a linear combination of discretized differential operators for some f(x) that represents the values of a function f on discretized x. To solve the system under boundary conditions on v, one can construct and solve the following extended system:beginbmatrix\noverlineL \nB\nendbmatrix \nv(overlinex) = \nbeginbmatrix\nf(x) \nb\nendbmatrix with M_E by (M+2) matrix B and M_E-length vector b that represent the current boundary conditions, where M_E is the number of boundary conditions to be applied. For instance, to apply reflecting barrier conditions v(x_min) = v(x_max) = 0, one can useB = beginbmatrix\n1  -1  0  dots  0  0  0 \n0  0  0  dots  0  -1  1\nendbmatrix_2 times (M+2) quad \nb = beginbmatrix\n0 \n0\nendbmatrixLikewise, for mixed boundary conditions v(x_min) = underlinexi and v(x_max) = overlinexi, one can useB = beginbmatrix\n1  -1  0  dots  0  0  0 \n0  0  0  dots  0  -1  1\nendbmatrix_2 times (M+2) quad \nb = beginbmatrix\n underlinexi Delta \noverlinexi Delta\nendbmatrix"
},

{
    "location": "derivations/#Applications-1",
    "page": "Derivations",
    "title": "Applications",
    "category": "section",
    "text": ""
},

{
    "location": "derivations/#Hamilton–Jacobi–Bellman-equations-(HJBE)-1",
    "page": "Derivations",
    "title": "Hamilton–Jacobi–Bellman equations (HJBE)",
    "category": "section",
    "text": "Consider solving for v from the following optimal control problemv(x_0) = max_ alpha(t)  _t geq 0  int_0^infty e^-rho t r( x(t) alpha(t )) dtwith the law of motion for the state dx = mu dt + sigma dW for some constant mu geq 0 and sigma geq 0 with x(0) = x_0.Let alpha^*(t) be the optimal solution. Suppose that r under alpha^*(t) can be expressed in terms of state variables, r^* (x). Then, the HJBE yields\\rho v(x) = r^*(x) +  \\mu  \\partial_{x} v(x) + \\dfrac{\\sigma^2}{2} \\partial_{xx} v(x)In terms of differential operators, one can rewrite the equation as(\\rho - \\tilde{L}) v(x) = r^*(x)where \\tilde{L} = \\mu \\partial_{x} + (\\sigma^2/2) \\partial_{xx}By descretizing the space of x, one can solve the corresponding system by using discretized operators for partial_x (L_1+), partial_xx (L_2) on some grids of length M, x_i_i=1^M:L = mu L_1+ + dfracsigma^22 L_2so that v under the optimal plan can be computed by solving the following discretized system of equations:(rho I - L) v = r^*where v and r^* are M-vectors whose ith elements are v(x_i) and r^*(x_i), respectively."
},

{
    "location": "derivations/#Kolmogorov-forward-equations-(KFE)-under-diffusion-process-1",
    "page": "Derivations",
    "title": "Kolmogorov forward equations (KFE) under diffusion process",
    "category": "section",
    "text": "Let g(x t) be the distribution of x at time t from the example above. By the Kolmogorov forward equation, the following PDE holds:\\partial_{t} g(x, t) = - \\mu \\partial_{x}  g(x,t) + \\dfrac{\\sigma^2}{2} \\partial_{xx} g(x,t)"
},

{
    "location": "derivations/#Stationary-distributions-1",
    "page": "Derivations",
    "title": "Stationary distributions",
    "category": "section",
    "text": "The stationary distribution g^*(x) satisfies0 = - mu partial_x g^*(x) + dfracsigma^22 partial_xx g^*(x)which can be rewritten as tildeL^* g(x) = 0where tildeL^* =  - mu partial_x + (sigma^22) partial_xxBy descretizing the space of x, one can solve the corresponding system by using discretized operators for tildeL^*. Note that the operator for the KFE in the original equation is the adjoint operator of the operator for the HJBE, tildeL, and the correct discretization scheme for L^* is, analogously, done by taking the transpose of the discretized operator for HJBE, L (See Gabaix et al., 2016). Hence, one can find the stationary distribution by solving the following discretized system of equations:L^T g = 0 where L^T is the transpose of L and g is an M-vector whose element is g(x_i) such that sum_i=1^M g(x_i) = 1."
},

{
    "location": "derivations/#Full-dynamics-of-distributions-1",
    "page": "Derivations",
    "title": "Full dynamics of distributions",
    "category": "section",
    "text": "One can also solve the full PDE in KFE equation, given an initial distribution g(x 0). After discretization, note that \\eqref{eq:kfe} can be rewritten asdotg(t) = L^T g(t)where dotg(t) is an M-vector whose ith element is partial_t g(x_i t), which can be efficently solved by a number of differential equation solvers available in public, including DifferentialEquations.jl."
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{Mixed,Mixed},DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x) by length(x) matrix under mixed boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   0.0   ⋅\n -1.0   1.0  0.0\n   ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅\n  0.0  -1.0   1.0\n   ⋅    0.0  -1.0\n\njulia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n 0.0   1.0    ⋅\n 1.0  -2.0   1.0\n  ⋅    1.0  -2.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{Reflecting,Reflecting},DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x) by length(x) matrix under reflecting boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> DifferentialOperator(x, (Reflecting(), Reflecting()), BackwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0   ⋅\n -1.0   1.0  0.0\n   ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x, (Reflecting(), Reflecting()), ForwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0   ⋅\n  0.0  -1.0  1.0\n   ⋅    0.0  0.0\n\njulia> DifferentialOperator(x, (Reflecting(), Reflecting()), CentralSecondDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅\n  1.0  -2.0   1.0\n   ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.ExtensionDifferentialOperator-Tuple{Any,DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.ExtensionDifferentialOperator",
    "category": "method",
    "text": "ExtensionDifferentialOperator(x, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x) by length(x) + 2 matrix whose first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively under no boundary condition using finite difference method specified by method.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> ExtensionDifferentialOperator(x, BackwardFirstDifference())\n3×5 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [1, 3]  =  0.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [2, 4]  =  0.0\n  [3, 4]  =  1.0\n\njulia> ExtensionDifferentialOperator(x, ForwardFirstDifference())\n3×5 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 2]  =  -1.0\n  [2, 2]  =  0.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -1.0\n  [3, 3]  =  0.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -1.0\n  [3, 5]  =  1.0\n\njulia> ExtensionDifferentialOperator(x, CentralSecondDifference())\n3×5 SparseArrays.SparseMatrixCSC{Float64,Int64} with 9 stored entries:\n  [1, 1]  =  1.0\n  [1, 2]  =  -2.0\n  [2, 2]  =  1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -2.0\n  [3, 3]  =  1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -2.0\n  [3, 5]  =  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₁₊-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₁₊",
    "category": "method",
    "text": "L̄₁₊(x)\n\nReturns a discretized first-order differential operator of length(x) by length(x) + 2 matrix using forward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> Array(L̄₁₊(x))\n3×5 Array{Float64,2}:\n 0.0  -1.0   1.0   0.0  0.0\n 0.0   0.0  -1.0   1.0  0.0\n 0.0   0.0   0.0  -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₁₋-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₁₋",
    "category": "method",
    "text": "L̄₁₋(x)\n\nReturns a discretized first-order differential operator of length(x) by length(x) + 2 matrix using backward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> Array(L̄₁₋(x))\n3×5 Array{Float64,2}:\n -1.0   1.0   0.0  0.0  0.0\n  0.0  -1.0   1.0  0.0  0.0\n  0.0   0.0  -1.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₂-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₂",
    "category": "method",
    "text": "L̄₂(x)\n\nReturns a discretized second-order differential operator of length(x) by length(x) + 2 matrix using central difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> Array(L̄₂(x))\n3×5 Array{Float64,2}:\n 1.0  -2.0   1.0   0.0  0.0\n 0.0   1.0  -2.0   1.0  0.0\n 0.0   0.0   1.0  -2.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊",
    "category": "method",
    "text": "L₁₊(x, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of length(x) by length(x) matrix using forward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L₁₊(x, (Reflecting(), Reflecting()))\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0   ⋅\n  0.0  -1.0  1.0\n   ⋅    0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋",
    "category": "method",
    "text": "L₁₋(x, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of length(x) by length(x) matrix using backward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L₁₋(x, (Reflecting(), Reflecting()))\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0   ⋅\n -1.0   1.0  0.0\n   ⋅   -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂",
    "category": "method",
    "text": "L₂(x, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized second-order differential operator of length(x) by length(x) matrix using central difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L₂(x, (Reflecting(), Reflecting()))\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅\n  1.0  -2.0   1.0\n   ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.x̄-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.x̄",
    "category": "method",
    "text": "x̄(x)\n\nReturns an extended grid of length length(x)+2 given grid x.\n\nThe first and last elements of the returned extended grid represent the ghost nodes just before x[1] and x[end] respectively.\n\njulia> x = 1:3\n1:3\n\njulia> x̄(x)\n5-element Array{Int64,1}:\n 0\n 1\n 2\n 3\n 4\n\njulia> x = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> x̄(x)\n5-element Array{Float64,1}:\n 0.5\n 1.0\n 1.5\n 1.7\n 1.9\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}
