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
    "text": "Consider solving for v from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):rho v(x) = f(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v under the reflecting barrier conditions v(0) = v(1) = 0 on M-size discretized grids, one can run the following code:# import LinearAlgebra package (for diagonal and identity matrices)\nusing LinearAlgebra \n# setup \nf(x) = x^2 \nμ = -0.1 # constant negative drift\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid\nx = range(0.0, 1.0, length = (M+2)) # grid\n\n# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nL = I * ρ - μ*L₁₋(x, bc) - σ^2 / 2 * L₂(̄, bc)\n## solve the value function\nv = L \\ f.(x) Note that the code above uses differential operators with reflecting boundary conditions applied.  One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute v:# import SparseArrays package (for identity matrix and spzeros)\nusing SparseArrays\n\n# differential operators on extended nodes\nL̄ₓ = μ*L̄₁₋(x) + σ^2 / 2 * L̄₂(x)\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])\nb = [0.0; 0.0] \n\n# form bellman equation on extension\nL̄ = [spzeros(M) ρ*I spzeros(M)] - L̄ₓ\n\n# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L̄; B] \\ [f.(x); b]\n\n# extract the interior (is identical with `v` above)\nv =  v̄[2:end-1] Here is a plot for v:using Plots\nplot(x, v, lw = 4, label = \"v\")(Image: plot-hjbe-both-reflecting)"
},

{
    "location": "#Solving-HJBE-with-absorbing-barrier-conditions-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with absorbing barrier conditions",
    "category": "section",
    "text": "Instead of having the reflecting barrier conditions on both lower bound and upper bound v(0) = v(1) = 0 as above, one can impose an absorbing barrier condition as well. To solve v under the reflecting barrier conditions v(0) = S (absorbing barrier on lower bound) for some S and v(1) = 0 (reflecting barrier on upper bound), one can construct B and b for the boundary conditions as follows:# define S\nS = 3.0 \n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[0; 1; zeros(M)] [zeros(M); -1; 1]])\nb = [S; 0.0];and solve v:# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L̄; B] \\ [f.(x); b]\n\n# extract the interior (is identical with `v` above)\nv =  v̄[2:end-1] Note that this can be alternatively done by Here is a plot for v:plot(x, v, lw = 4, label = \"v\")(Image: plot-hjbe-lb-absorbing-ub-reflecting)"
},

{
    "location": "#Solving-HJBE-with-state-dependent-drifts-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with state-dependent drifts",
    "category": "section",
    "text": "One can also deploy upwind schemes when drift variable is not constant. Consider solving for v from the following Bellman equation:rho v(x) = f(x) + mu(x) partial_x v(x) + fracsigma^22 partial_xx v(x)associated with the diffusion processdx = mu(x) dt + sigma dWfor some constant rho sigma  0 and mu(x) = -x. Note that mu(x) depends on states. The following code will solve v using upwind schemes, with the reflecting barrier conditions v(0) = v(1) = 0 applied:# setup \nf(x) = x^2 \nμ(x) = -x # drift depends on state\nσ = 1.0\nρ = 0.05\nM = 100 # size of grid\nx = range(-1.0, 1.0, length = (M+2))\n\nbc = (Reflecting(), Reflecting())\n\n# Define first order differential operator using upwind scheme\nL₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋(x, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊(x, bc)\n\n# Define linear operator using upwind schemes\nL = L₁ - σ^2 / 2 * L₂(x,bc)\n\n# solve the value function\nv = (I * ρ - L) \\ f.(x) "
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
    "text": "DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) matrix under mixed boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n Inf     0.0    ⋅    ⋅\n  -1.0   1.0   0.0   ⋅\n    ⋅   -1.0   1.0  0.0\n    ⋅     ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  0.0  -1.0   1.0    ⋅\n   ⋅    0.0  -1.0   1.0\n   ⋅     ⋅    0.0  -0.5\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -Inf     1.0    ⋅     ⋅\n    1.0  -2.0   1.0    ⋅\n     ⋅    1.0  -2.0   1.0\n     ⋅     ⋅    1.0  -1.5\n\n\n\n\n\n"
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
    "location": "api/#SimpleDifferentialOperators.interiornodes-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.interiornodes",
    "category": "method",
    "text": "interiornodes(x̄, bc)\n\nReturns an interior grid corresponding to the boundary condition bc given extended grid x̄.\n\njulia> x̄ = 0:5\n0:5\n\njulia> interiornodes(x̄, (Reflecting(), Reflecting()))\n1:4\n\njulia> x̄ = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> interiornodes(x̄, (Mixed(1.0), Mixed(1.0)))\n1-element Array{Float64,1}:\n 1.5\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.interiornodes-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.interiornodes",
    "category": "method",
    "text": "interiornodes(x̄)\n\nReturns an interior grid of length length(x̄)-2 given extended grid x̄.\n\njulia> x̄ = 0:5\n0:5\n\njulia> interiornodes(x̄)\n1:4\n\njulia> x̄ = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> interiornodes(x̄)\n1-element Array{Float64,1}:\n 1.5\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}