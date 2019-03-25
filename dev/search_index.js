var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "﻿SimpleDifferentialOperators.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#[SimpleDifferentialOperators.jl](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/)-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "﻿SimpleDifferentialOperators.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Installation-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Installation",
    "category": "section",
    "text": "To install, run] add SimpleDifferentialOperatorsNote that this requires Julia 1.1 or later."
},

{
    "location": "#Usage-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Usage",
    "category": "section",
    "text": ""
},

{
    "location": "#Solving-HJBE-with-constant-drifts-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with constant drifts",
    "category": "section",
    "text": "Consider solving for v from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):rho v(x) = f(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v under the reflecting barrier conditions v(0) = v(1) = 0 on M-size discretized grids, one can run the following code:# import LinearAlgebra package (for diagonal and identity matrices)\nusing LinearAlgebra \n# setup \nf(x) = x^2 \nμ = -0.1 # constant negative drift\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid\nx̄ = range(0.0, 1.0, length = (M+2)) # grid\nx = interior(x̄) # interior nodes\n# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nL = I * ρ - μ*L₁₋(x̄, bc) - σ^2 / 2 * L₂(x̄, bc)\n## solve the value function\nv = L \\ f.(x) Note that the code above uses differential operators with reflecting boundary conditions applied.  One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute v:# import SparseArrays package (for identity matrix and spzeros)\nusing SparseArrays\n\n# differential operators on extended nodes\nL̄ₓ = μ*L̄₁₋(x̄) + σ^2 / 2 * L̄₂(x̄)\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])\nb = [0.0; 0.0] \n\n# form bellman equation on extension\nL̄ = [spzeros(M) ρ*I spzeros(M)] - L̄ₓ\n\n# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L̄; B] \\ [f.(x); b]\n\n# extract the interior (is identical with `v` above)\nv =  v̄[2:end-1] Here is a plot for v:using Plots\nplot(x, v, lw = 4, label = \"v\")(Image: plot-hjbe-both-reflecting)"
},

{
    "location": "#Solving-HJBE-with-absorbing-barrier-conditions-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with absorbing barrier conditions",
    "category": "section",
    "text": "Instead of having the reflecting barrier conditions on both lower bound and upper bound v(0) = v(1) = 0 as above, one can impose an absorbing barrier condition as well. To solve v under the reflecting barrier conditions v(0) = S (absorbing barrier on lower bound) for some S and v(1) = 0 (reflecting barrier on upper bound), one can construct B and b for the boundary conditions as follows:# define S\nS = 0.0 \n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[0; 1; zeros(M)] [zeros(M); -1; 1]])\nb = [S; 0.0];and solve v:# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L̄; B] \\ [f.(x); b]\n\n# extract the interior (is identical with `v` above)\nv =  v̄[2:end-1] Note that this can be alternatively done by Here is a plot for v:plot(x, v, lw = 4, label = \"v\")(Image: plot-hjbe-lb-absorbing-ub-reflecting)"
},

{
    "location": "#Solving-HJBE-with-state-dependent-drifts-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with state-dependent drifts",
    "category": "section",
    "text": "One can also deploy upwind schemes when drift variable is not constant. Consider solving for v from the following Bellman equation:rho v(x) = f(x) + mu(x) partial_x v(x) + fracsigma^22 partial_xx v(x)associated with the diffusion processdx = mu(x) dt + sigma dWfor some constant rho sigma  0 and mu(x) = -x. Note that mu(x) depends on states. The following code will solve v using upwind schemes, with the reflecting barrier conditions v(0) = v(1) = 0 applied:# setup \nf(x) = x^2 \nμ(x) = -x # drift depends on state\nσ = 1.0\nρ = 0.05\nM = 100 # size of grid\nx̄ = range(-1.0, 1.0, length = (M+2))\nx = interior(x̄) # interior nodes\n\nbc = (Reflecting(), Reflecting())\n\n# Define first order differential operator using upwind scheme\nL₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋(x̄, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊(x̄, bc)\n\n# Define linear operator using upwind schemes\nL = L₁ - σ^2 / 2 * L₂(x̄,bc)\n\n# solve the value function\nv = (I * ρ - L) \\ f.(x) "
},

{
    "location": "#Finding-stationary-distribution-from-the-Kolmogorov-forward-equation-(KFE)-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Finding stationary distribution from the Kolmogorov forward equation (KFE)",
    "category": "section",
    "text": "One can also compute the stationary distribution of the state x above from the corresponding KFE:partial_t g(xt) = - mu(x) partial_x g(x t) + fracsigma^22 partial_xx g(xt)by taking partial_t g(xt) = 0, i.e., solving g from the L^* g(x) = 0 whereL^* = - mu(x) partial_x + fracsigma^22 partial_xxBy descretizing the space of x, one can solve the corresponding system by using discretized operators for L^*. Note that the operator for the KFE in the original equation is the adjoint operator of the operator for the HJBE, L, and the correct discretization scheme for L^* is, analogously, done by taking the transpose of the discretized operator for HJBE, L (See Gabaix et al., 2016 and Achdou et al., 2017). Hence, one can find the stationary distribution by solving the following discretized system of equations:L^T g = 0such that the sum of g is one. This can be found by finding a non-trivial eigenvector for L^T  associated with the eigenvalue of zero:using Arpack # library for extracting eigenvalues and eigenvectors\n\n# extract eigenvalues and eigenvectors, smallest eigenval in magintute first\nλ, ϕ = eigs(transpose(L), which = :SM); \n# extract the very first eigenvector (associated with the smallest eigenvalue)\ng_ss = real.(ϕ[:,1]);\n# normalize it\ng_ss = g_ss / sum(g_ss)Using L from the state-dependent drift example above, this results in the following stationary distribution:plot(x, g_ss, lw = 4, label = \"g_ss\")(Image: plot-stationary-dist)"
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
    "text": "We prepared a gallery of notebooks to showcase the applicability of the package to different economic problems.Optimal Stopping Problems. You can download the notebook here.\nComputational Appendix: Optimal Stopping Problems. You can download the notebook here. This notebook presents solution to the problem in Julia, using a number of solvers and modeling tools."
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
    "text": "Let overlinex = x_i_i=0^M+1 be an extended grid of discretized M-length of grids on x with end points x_0 = x_min and x_M+1 = x_max, and x = x_i_i=1^M be the corresponding interior grid. Also, throughout the section, we consider regular grids, i.e., x_i+1 - x_i = Delta for some constant Delta  0 for all i = 0  M. Given a real-valued function v, let v(x) be the M-length vector whose ith element is v(x_i). The goal is to construct a matrix L such that L v(x) represents the first-order or second-order derivative of v on x under some boundary conditions."
},

{
    "location": "derivations/#Mixed-Boundary-Values-1",
    "page": "Derivations",
    "title": "Mixed Boundary Values",
    "category": "section",
    "text": "Under M-length of grids on v with end points x_min  x_max with mixed boundary value conditions ofbeginalign\nunderlinexi v(x_min) + nabla v(x_min) = 0 \noverlinexi v(x_max) + nabla v(x_max) = 0\nendalignNote that this can be extended to reflecting boundary conditions by assigning underlinexi = 0 and overlinexi = 0. We use the following discretization schemes:L_1- equiv frac1Deltabeginpmatrix\n1 - (1 + underlinexi Delta) 00dots000\n-110dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Mtimes ML_1+ equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots0-11\n000cdots00-1+(1-overlinexi Delta)\nendpmatrix_Mtimes MlabeleqL-1-plus-regular L_2 equiv frac1Delta^2beginpmatrix\n-2 + (1 + underlinexi Delta) 10dots000\n1-21dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots1-21\n000cdots01-2 + (1- overlinexi Delta)\nendpmatrix_Mtimes MlabeleqL-2-regularwhich represent the backward first order, foward first order, and central second order differential operators respectively."
},

{
    "location": "derivations/#Applying-boundary-conditions-with-operators-on-extended-grids-1",
    "page": "Derivations",
    "title": "Applying boundary conditions with operators on extended grids",
    "category": "section",
    "text": "Boundary conditions can be applied manually by using operators on extended grids. Define v(overlinex) as (M+2)-vector whose ith element is overline x_i-1. We can then define the following operators on overlinex:overlineL_1- equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots100\n000cdots-110\nendpmatrix_Mtimes (M+2)overlineL_1+ equiv frac1Deltabeginpmatrix\n0-11dots000\n00-1dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Mtimes (M+2) overlineL_2 equiv frac1Delta^2beginpmatrix\n-12-1dots000\n0-12dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots2-10\n000cdots-121\nendpmatrix_Mtimes (M+2)Suppose that we want to solve a system L v(x) = f(x)  where L is a linear combination of discretized differential operators for some f(x) that represents the values of a function f on discretized x. To solve the system under boundary conditions on v, one can construct and solve the following extended system:beginbmatrix\noverlineL \nB\nendbmatrix \nv(overlinex) = \nbeginbmatrix\nf(x) \nb\nendbmatrix with M_E by (M+2) matrix B and M_E-length vector b that represent the current boundary conditions, where M_E is the number of boundary conditions to be applied. For instance, to apply reflecting barrier conditions v(x_min) = v(x_max) = 0, one can useB = beginbmatrix\n-1  1  0  dots  0  0  0 \n0  0  0  dots  0  -1  1\nendbmatrix_2 times (M+2) quad \nb = beginbmatrix\n0 \n0\nendbmatrixLikewise, for non-homogenous boundary conditions v(x_min) = underlinexi and v(x_max) = overlinexi, one can useB = beginbmatrix\n-1  1  0  dots  0  0  0 \n0  0  0  dots  0  -1  1\nendbmatrix_2 times (M+2) quad \nb = beginbmatrix\n underlinexi Delta \noverlinexi Delta\nendbmatrix"
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
    "text": "Consider solving for v from the following optimal control problemv(x_0) = max_ alpha(t)  _t geq 0  int_0^infty e^-rho t r( x(t) alpha(t )) dtwith the law of motion for the state dx = mu dt + sigma dW for some constant mu geq 0 and sigma geq 0 with x(0) = x_0.Let alpha^*(t) be the optimal solution. Suppose that r under alpha^*(t) can be expressed in terms of state variables, r^* (x). Then, the HJBE yieldsrho v(x) = r^*(x) +  mu  partial_x v(x) + dfracsigma^22 partial_xx v(x)In terms of differential operators, one can rewrite the equation as(rho - tildeL) v(x) = r^*(x)where tildeL = mu partial_x + (sigma^22) partial_xxBy descretizing the space of x, one can solve the corresponding system by using discretized operators for partial_x (L_1+), partial_xx (L_2) on some grids of length M, x_i_i=1^M:L = mu L_1+ + dfracsigma^22 L_2so that v under the optimal plan can be computed by solving the following discretized system of equations:(rho I - L) v = r^*where v and r^* are M-vectors whose ith elements are v(x_i) and r^*(x_i), respectively."
},

{
    "location": "derivations/#Kolmogorov-forward-equations-(KFE)-under-diffusion-process-1",
    "page": "Derivations",
    "title": "Kolmogorov forward equations (KFE) under diffusion process",
    "category": "section",
    "text": "Let g(x t) be the distribution of x at time t from the example above. By the Kolmogorov forward equation, the following PDE holds:partial_t g(x t) = - mu partial_x  g(xt) + dfracsigma^22 partial_xx g(xt)"
},

{
    "location": "derivations/#Stationary-distributions-1",
    "page": "Derivations",
    "title": "Stationary distributions",
    "category": "section",
    "text": "The stationary distribution g^*(x) satisfies0 = - mu partial_x g^*(x) + dfracsigma^22 partial_xx g^*(x)which can be rewritten as tildeL^* g(x) = 0where tildeL^* =  - mu partial_x + (sigma^22) partial_xxBy descretizing the space of x, one can solve the corresponding system by using discretized operators for tildeL^*. Note that the operator for the KFE in the original equation is the adjoint operator of the operator for the HJBE, tildeL, and the correct discretization scheme for tildeL^* is, analogously, done by taking the transpose of the discretized operator for HJBE, L (See Gabaix et al., 2016 and Achdou et al., 2017). Hence, one can find the stationary distribution by solving the following discretized system of equations:L^T g = 0 where L^T is the transpose of L and g is an M-vector whose element is g(x_i) such that sum_i=1^M g(x_i) = 1."
},

{
    "location": "derivations/#Full-dynamics-of-distributions-1",
    "page": "Derivations",
    "title": "Full dynamics of distributions",
    "category": "section",
    "text": "One can also solve the full PDE in KFE equation, given an initial distribution g(x 0). After discretization, note that the KFE can be rewritten asdotg(t) = L^T g(t)where dotg(t) is an M-vector whose ith element is partial_t g(x_i t), which can be efficently solved by a number of differential equation solvers available in public, including DifferentialEquations.jl."
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
    "text": "DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) matrix under mixed boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   0.0    ⋅    ⋅\n -1.0   1.0   0.0   ⋅\n   ⋅   -1.0   1.0  0.0\n   ⋅     ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  0.0  -1.0   1.0    ⋅\n   ⋅    0.0  -1.0   1.0\n   ⋅     ⋅    0.0  -1.0\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n 0.0   1.0    ⋅     ⋅\n 1.0  -2.0   1.0    ⋅\n  ⋅    1.0  -2.0   1.0\n  ⋅     ⋅    1.0  -2.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{Reflecting,Reflecting},DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "DifferentialOperator(x̄, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) matrix under reflecting boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), BackwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0    ⋅    ⋅\n -1.0   1.0   0.0   ⋅\n   ⋅   -1.0   1.0  0.0\n   ⋅     ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), ForwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅    ⋅\n  0.0  -1.0   1.0   ⋅\n   ⋅    0.0  -1.0  1.0\n   ⋅     ⋅    0.0  0.0\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), CentralSecondDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  1.0  -2.0   1.0    ⋅\n   ⋅    1.0  -2.0   1.0\n   ⋅     ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.ExtensionDifferentialOperator-Tuple{Any,DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.ExtensionDifferentialOperator",
    "category": "method",
    "text": "ExtensionDifferentialOperator(x̄, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) + 2 matrix whose first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively under no boundary condition using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> ExtensionDifferentialOperator(x̄, BackwardFirstDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 11 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [1, 3]  =  0.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [2, 4]  =  0.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -1.0\n  [3, 5]  =  0.0\n  [4, 5]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, ForwardFirstDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 11 stored entries:\n  [1, 2]  =  -1.0\n  [2, 2]  =  0.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -1.0\n  [3, 3]  =  0.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -1.0\n  [4, 4]  =  0.0\n  [3, 5]  =  1.0\n  [4, 5]  =  -1.0\n  [4, 6]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, CentralSecondDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 12 stored entries:\n  [1, 1]  =  1.0\n  [1, 2]  =  -2.0\n  [2, 2]  =  1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -2.0\n  [3, 3]  =  1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -2.0\n  [4, 4]  =  1.0\n  [3, 5]  =  1.0\n  [4, 5]  =  -2.0\n  [4, 6]  =  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₁₊-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₁₊",
    "category": "method",
    "text": "L̄₁₊(x̄)\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) + 2 matrix using forward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> Array(L̄₁₊(x̄))\n4×6 Array{Float64,2}:\n 0.0  -1.0   1.0   0.0   0.0  0.0\n 0.0   0.0  -1.0   1.0   0.0  0.0\n 0.0   0.0   0.0  -1.0   1.0  0.0\n 0.0   0.0   0.0   0.0  -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₁₋-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₁₋",
    "category": "method",
    "text": "L̄₁₋(x̄)\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) + 2 matrix using backward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively.\n\nExamples\n\njulia> x̄ = 1:3\n1:3\n\njulia> Array(L̄₁₋(x̄))\n1×3 Array{Float64,2}:\n -1.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₂-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₂",
    "category": "method",
    "text": "L̄₂(x̄)\n\nReturns a discretized second-order differential operator of length(x̄) by length(x̄) + 2 matrix using central difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively.\n\nExamples\n\njulia> x̄ = 0:5\n0:5 \n\njulia> Array(L̄₂(x̄))\n4×6 Array{Float64,2}:\n 1.0  -2.0   1.0   0.0   0.0  0.0\n 0.0   1.0  -2.0   1.0   0.0  0.0\n 0.0   0.0   1.0  -2.0   1.0  0.0\n 0.0   0.0   0.0   1.0  -2.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊",
    "category": "method",
    "text": "L₁₊(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) matrix using forward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₊(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅    ⋅\n  0.0  -1.0   1.0   ⋅\n   ⋅    0.0  -1.0  1.0\n   ⋅     ⋅    0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋",
    "category": "method",
    "text": "L₁₋(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) matrix using backward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₋(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0    ⋅    ⋅\n -1.0   1.0   0.0   ⋅\n   ⋅   -1.0   1.0  0.0\n   ⋅     ⋅   -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂",
    "category": "method",
    "text": "L₂(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized second-order differential operator of length(x̄) by length(x̄) matrix using central difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₂(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  1.0  -2.0   1.0    ⋅\n   ⋅    1.0  -2.0   1.0\n   ⋅     ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.interior-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.interior",
    "category": "method",
    "text": "interior(x̄, bc)\n\nReturns an interior grid corresponding to the boundary condition bc given extended grid x̄.\n\njulia> x̄ = 0:5\n0:5\n\njulia> interior(x̄, (Reflecting(), Reflecting()))\n1:4\n\njulia> x̄ = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> interior(x̄, (Mixed(1.0), Mixed(1.0)))\n1-element Array{Float64,1}:\n 1.5\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.interior-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.interior",
    "category": "method",
    "text": "interior(x̄)\n\nReturns an interior grid of length length(x̄)-2 given extended grid x̄.\n\njulia> x̄ = 0:5\n0:5\n\njulia> interior(x̄)\n1:4\n\njulia> x̄ = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> interior(x̄)\n1-element Array{Float64,1}:\n 1.5\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}
