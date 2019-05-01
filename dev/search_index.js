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
    "text": "Detailed derivations and more applications can be found here."
},

{
    "location": "#Solving-HJBE-with-constant-drifts-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with constant drifts",
    "category": "section",
    "text": "Consider solving for v from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):rho v(x) = f(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v under the reflecting barrier conditions v(0) = v(1) = 0 on M-size discretized grids, one can run the following code:using LinearAlgebra, SimpleDifferentialOperators\n# setup\nf(x) = x^2\nμ = -0.1 # constant negative drift\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid (interior points)\n\nx̄ = range(0.0, 1.0, length = (M+2))\nx = interiornodes(x̄) # i.e., x̄[2:end-1]\n\n# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nL_bc = I * ρ - μ*L₁₋bc(x̄, bc) - σ^2 / 2 * L₂bc(x̄, bc)\n\n# solve the value function\nv = L_bc \\ f.(x)Note that the interior solution v does not the values of v at the boundary, i.e., v(0) and v(1). To extend the interior solution to the boundary points, one can call extrapolatetoboundary as follows:̄v = extrapolatetoboundary(x̄, v, bc);Here is a complete plot for v:using Plots\nplot(x̄, v̄, lw = 4, label = \"v\")(Image: plot-hjbe-both-reflecting)Note that the code above uses differential operators on the interior nodes with reflecting boundary conditions applied. One can alternatively use operators on extended nodes (extended operators) and stack them with matrices for boundary conditions to compute v:# import SparseArrays package (for identity matrix and spzeros)\nusing SparseArrays\n\n# differential operators on extended nodes\nLₓ = μ*L₁₋(x̄) + σ^2 / 2 * L₂(x̄)\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])\nb = [0.0; 0.0]\n\n# form bellman equation on extension\nL = [spzeros(M) ρ*I spzeros(M)] - Lₓ\n\n# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L; B] \\ [f.(x); b]\n\n# extract the interior (is identical with `v` above)\nv =  v̄[2:end-1]"
},

{
    "location": "#Solving-HJBE-with-absorbing-barrier-conditions-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with absorbing barrier conditions",
    "category": "section",
    "text": "Instead of having the reflecting barrier conditions on both lower bound and upper bound v(0) = v(1) = 0 as above, one can impose an absorbing barrier condition as well. To solve v under the reflecting barrier conditions v(0) = S (absorbing barrier on lower bound) for some S and v(1) = 0 (reflecting barrier on upper bound), one can construct B and b for the boundary conditions as follows:# define S\nS = 3.0\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])\nb = [S; 0.0];and solve v:# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L; B] \\ [f.(x); b]Note that this can be alternatively done byHere is a plot for v:plot(x̄, v̄, lw = 4, label = \"v\")(Image: plot-hjbe-lb-absorbing-ub-reflecting)"
},

{
    "location": "#Solving-HJBE-with-state-dependent-drifts-1",
    "page": "﻿SimpleDifferentialOperators.jl",
    "title": "Solving HJBE with state-dependent drifts",
    "category": "section",
    "text": "One can also deploy upwind schemes when drift variable is not constant. Consider solving for v from the following Bellman equation:rho v(x) = f(x) + mu(x) partial_x v(x) + fracsigma^22 partial_xx v(x)associated with the diffusion processdx = mu(x) dt + sigma dWfor some constant rho sigma  0 and mu(x) = -x. Note that mu(x) depends on states. The following code will solve v using upwind schemes, with the reflecting barrier conditions v(0) = v(1) = 0 applied:# setup\nf(x) = x^2\nμ(x) = -x # drift depends on state\nσ = 1.0\nρ = 0.05\nM = 100 # size of grid\n\nx̄ = range(-1., 1., length = M + 2)\nx = interiornodes(x̄) # i.e., x̄[2:end-1]\n\nbc = (Reflecting(), Reflecting())\n\n# Define first order differential operator using upwind scheme\nL₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋bc(x̄, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊bc(x̄, bc)\n\n# Define linear operator using upwind schemes\nL_x = L₁ - σ^2 / 2 * L₂bc(x̄, bc)\nL_bc = I * ρ - L_x\n\n# solve the value function\nv = L_bc \\ f.(x)"
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
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{Mixed,Mixed},BackwardFirstDifference}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "DifferentialOperator(x̄, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) matrix under mixed boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 7 stored entries:\n  [1, 1]  =  Inf\n  [2, 1]  =  -1.0\n  [2, 2]  =  1.0\n  [3, 2]  =  -1.0\n  [3, 3]  =  1.0\n  [4, 3]  =  -1.0\n  [4, 4]  =  1.0\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 7 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -0.5\n\njulia> DifferentialOperator(x̄, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 10 stored entries:\n  [1, 1]  =  -Inf\n  [2, 1]  =  1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -2.0\n  [3, 2]  =  1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -2.0\n  [4, 3]  =  1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -1.5\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{Reflecting,Reflecting},Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "DifferentialOperator(x̄, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) matrix under reflecting boundary conditions from bc using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), BackwardFirstDifference())\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 7 stored entries:\n  [1, 1]  =  0.0\n  [2, 1]  =  -1.0\n  [2, 2]  =  1.0\n  [3, 2]  =  -1.0\n  [3, 3]  =  1.0\n  [4, 3]  =  -1.0\n  [4, 4]  =  1.0\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), ForwardFirstDifference())\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 7 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  0.0\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), CentralSecondDifference())\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 10 stored entries:\n  [1, 1]  =  -1.0\n  [2, 1]  =  1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -2.0\n  [3, 2]  =  1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -2.0\n  [4, 3]  =  1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.ExtensionDifferentialOperator-Tuple{AbstractRange,BackwardFirstDifference}",
    "page": "API",
    "title": "SimpleDifferentialOperators.ExtensionDifferentialOperator",
    "category": "method",
    "text": "ExtensionDifferentialOperator(x̄, method::DifferenceMethod)\n\nReturns a discretized differential operator of length(x̄) by length(x̄) + 2 matrix whose first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively under no boundary condition using finite difference method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> ExtensionDifferentialOperator(x̄, BackwardFirstDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -1.0\n  [4, 5]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, ForwardFirstDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 2]  =  -1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -1.0\n  [3, 5]  =  1.0\n  [4, 5]  =  -1.0\n  [4, 6]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, CentralSecondDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 12 stored entries:\n  [1, 1]  =  1.0\n  [1, 2]  =  -2.0\n  [2, 2]  =  1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -2.0\n  [3, 3]  =  1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -2.0\n  [4, 4]  =  1.0\n  [3, 5]  =  1.0\n  [4, 5]  =  -2.0\n  [4, 6]  =  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊",
    "category": "method",
    "text": "L₁₊(x̄)\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) + 2 matrix using forward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> Array(L₁₊(x̄))\n4×6 Array{Float64,2}:\n 0.0  -1.0   1.0   0.0   0.0  0.0\n 0.0   0.0  -1.0   1.0   0.0  0.0\n 0.0   0.0   0.0  -1.0   1.0  0.0\n 0.0   0.0   0.0   0.0  -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊bc-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊bc",
    "category": "method",
    "text": "L₁₊bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) matrix using forward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₊bc(x̄, (Reflecting(), Reflecting()))\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 7 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋",
    "category": "method",
    "text": "L₁₋(x̄)\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) + 2 matrix using backward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively.\n\nExamples\n\njulia> x̄ = 1:3\n1:3\n\njulia> Array(L₁₋(x̄))\n1×3 Array{Float64,2}:\n -1.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋bc-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋bc",
    "category": "method",
    "text": "L₁₋bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of length(x̄) by length(x̄) matrix using backward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₋bc(x̄, (Reflecting(), Reflecting()))\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 7 stored entries:\n  [1, 1]  =  0.0\n  [2, 1]  =  -1.0\n  [2, 2]  =  1.0\n  [3, 2]  =  -1.0\n  [3, 3]  =  1.0\n  [4, 3]  =  -1.0\n  [4, 4]  =  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂",
    "category": "method",
    "text": "L₂(x̄)\n\nReturns a discretized second-order differential operator of length(x̄) by length(x̄) + 2 matrix using central difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x̄[1] and x̄[end] respectively.\n\nExamples\n\njulia> x̄ = 0:5\n0:5 \n\njulia> Array(L₂(x̄))\n4×6 Array{Float64,2}:\n 1.0  -2.0   1.0   0.0   0.0  0.0\n 0.0   1.0  -2.0   1.0   0.0  0.0\n 0.0   0.0   1.0  -2.0   1.0  0.0\n 0.0   0.0   0.0   1.0  -2.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂bc-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂bc",
    "category": "method",
    "text": "L₂bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized second-order differential operator of length(x̄) by length(x̄) matrix using central difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₂bc(x̄, (Reflecting(), Reflecting()))\n4×4 SparseArrays.SparseMatrixCSC{Float64,Int64} with 10 stored entries:\n  [1, 1]  =  -1.0\n  [2, 1]  =  1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -2.0\n  [3, 2]  =  1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -2.0\n  [4, 3]  =  1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.extrapolatetoboundary-Tuple{Any,Any,Tuple{Mixed,Mixed}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.extrapolatetoboundary",
    "category": "method",
    "text": "extrapolatetoboundary(v, x̄, bc::Tuple{Mixed, Mixed})\n\nReturns a length(x̄)-vector whose 2:(length(x̄)-1) elements are v, the first and last element are extrapolated v on the boundaries of x̄ according to boundary conditions bc given.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> x = interiornodes(x̄)\n1:4\n\njulia> v = (x -> x^2).(x)\n4-element Array{Int64,1}:\n  1\n  4\n  9\n 16\n\njulia> extrapolatetoboundary(v, x̄, (Mixed(1), Mixed(1)))\n6-element Array{Float64,1}:\n Inf\n   1.0\n   4.0\n   9.0\n  16.0\n   8.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.extrapolatetoboundary-Tuple{Any,Any,Tuple{Reflecting,Reflecting}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.extrapolatetoboundary",
    "category": "method",
    "text": "extrapolatetoboundary(v, x̄, bc::Tuple{Reflecting, Reflecting})\n\nReturns a length(x̄)-vector whose 2:(length(x̄)-1) elements are v, the first and last element are extrapolated v on the boundaries of x̄ according to boundary conditions bc given.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> x = interiornodes(x̄)\n1:4\n\njulia> v = (x -> x^2).(x)\n4-element Array{Int64,1}:\n  1\n  4\n  9\n 16\n\njulia> extrapolatetoboundary(v, x̄, (Reflecting(), Reflecting()))\n6-element Array{Int64,1}:\n  1\n  1\n  4\n  9\n 16\n 16\n\n\n\n\n\n"
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
