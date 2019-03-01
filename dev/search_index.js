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
    "text": "Consider solving for v from the following equation:rho v(x) = f(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v on M-size discretized grids, one can run the following code:# import LinearAlgebra package (for identity matrices)\nusing LinearAlgebra \n# setup \nf(x) = x^2 \nμ = -0.1 # constant negative drift\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid\nx = range(0.0, 1.0, length = M) # grid\n\n# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nL = I * ρ - μ*L₁₋(x, bc) - σ^2 / 2 * L₂(x, bc)\n## solve the value function\nv_bc = L \\ f.(x) Note that the code above uses differential operators with reflecting boundary conditions applied.  One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute v:# operators without boundary conditions, adding extra two rows for boundary conditions\n## differential operators on extended nodes\nA = μ*L̄₁₋(x) + σ^2 / 2 * L̄₂(x)\n## matrix for boundary conditions\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]]) \n## stack them together\nL = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * 0.05 - A); B] \n## solve the value function with reflecting barrier bc (last two elements)\nv_bar = L \\ [f.(x); 0.0; 0.0] \n## extract the interior (is identical with `v_bc` above)\nv_interior = v_bar[2:end-1] "
},

{
    "location": "#Examples-1",
    "page": "﻿Installation",
    "title": "Examples",
    "category": "section",
    "text": "One can also deploy upwind schemes when drift variable is not constant. Consider solving for v from the following equation:rho v(x) = f(x) + mu(x) partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu(x) = -x. Note that mu(x) depends on states. The following code will solve v using upwind schemes:# setup \nf(x) = x^2 \nμ(x) = -x # drift depends on state\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid\nx = range(-1.0, 1.0, length = 100)\n\nbc = (Reflecting(), Reflecting())\n\n# Define first order differential operator using upwind scheme\nL₁ = min.(μ.(x), 0.0) .* L₁₋(x, bc) + max.(μ.(x), 0.0) .* L₁₊(x, bc)\n\n# Define linear operator using upwind schemes\nL = I * ρ - L₁ - σ^2 / 2 * L₂(x,bc)\n\n# solve the value function\nv_bc = L \\ f.(x) "
},

{
    "location": "formula/#",
    "page": "Formula",
    "title": "Formula",
    "category": "page",
    "text": ""
},

{
    "location": "formula/#Formula-1",
    "page": "Formula",
    "title": "Formula",
    "category": "section",
    "text": "Under M-length of grids on v with end points x_min  x_max with reflecting barrier conditions ofbeginalign\nunderlinexi v(x_min) + nabla v(x_min) = 0\noverlinexi v(x_max) + nabla v(x_max) = 0\nendalignWe use the following discretization schemes:L_1^- equiv frac1Deltabeginpmatrix\n1 - (1 + underlinexi Delta) 00dots000\n-110dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Ptimes PL_1^+ equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots0-11\n000cdots00-1+(1-overlinexi Delta)\nendpmatrix_Ptimes PlabeleqL-1-plus-regular L_2 equiv frac1Delta^2beginpmatrix\n-2 + (1 + underlinexi Delta) 10dots000\n1-21dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots1-21\n000cdots01-2 + (1- overlinexi Delta)\nendpmatrix_Ptimes PlabeleqL-2-regularwhich represent the backward first order, foward first order, and central second order differential operators respectively.Derivation, including formula for irregular grids, can be found here."
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
    "text": "`DifferentialOperator(x, bc::Tuple{Mixed, Mixed}, method::DifferenceMethod)`\n\nReturns a discretized differential operator of length(x) by length(x) matrix under mixed boundary conditions from bc using finite difference method specified by method. \n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), BackwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   0.0   ⋅\n -1.0   1.0  0.0\n   ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), ForwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅\n  0.0  -1.0   1.0\n   ⋅    0.0  -1.0\n\njulia> DifferentialOperator(x, (Mixed(1.0), Mixed(1.0)), CentralSecondDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n 0.0   1.0    ⋅\n 1.0  -2.0   1.0\n  ⋅    1.0  -2.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{Reflecting,Reflecting},DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "`DifferentialOperator(x, bc::Tuple{Reflecting, Reflecting}, method::DifferenceMethod)`\n\nReturns a discretized differential operator of length(x) by length(x) matrix under reflecting boundary conditions from bc using finite difference method specified by method. \n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> DifferentialOperator(x, (Reflecting(), Reflecting()), BackwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0   ⋅\n -1.0   1.0  0.0\n   ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x, (Reflecting(), Reflecting()), ForwardFirstDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0   ⋅\n  0.0  -1.0  1.0\n   ⋅    0.0  0.0\n\njulia> DifferentialOperator(x, (Reflecting(), Reflecting()), CentralSecondDifference())\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅\n  1.0  -2.0   1.0\n   ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.ExtensionDifferentialOperator-Tuple{Any,DifferenceMethod}",
    "page": "API",
    "title": "SimpleDifferentialOperators.ExtensionDifferentialOperator",
    "category": "method",
    "text": "`ExtensionDifferentialOperator(x, method::DifferenceMethod)`\n\nReturns a discretized differential operator of length(x) by length(x) + 2 matrix whose first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively under no boundary condition using finite difference method specified by method.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> ExtensionDifferentialOperator(x, BackwardFirstDifference())\n3×5 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [1, 3]  =  0.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [2, 4]  =  0.0\n  [3, 4]  =  1.0\n\njulia> ExtensionDifferentialOperator(x, ForwardFirstDifference())\n3×5 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 2]  =  -1.0\n  [2, 2]  =  0.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -1.0\n  [3, 3]  =  0.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -1.0\n  [3, 5]  =  1.0\n\njulia> ExtensionDifferentialOperator(x, CentralSecondDifference())\n3×5 SparseArrays.SparseMatrixCSC{Float64,Int64} with 9 stored entries:\n  [1, 1]  =  1.0\n  [1, 2]  =  -2.0\n  [2, 2]  =  1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -2.0\n  [3, 3]  =  1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -2.0\n  [3, 5]  =  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₁₊-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₁₊",
    "category": "method",
    "text": "`L̄₁₊(x)`\n\nReturns a discretized first-order differential operator of length(x) by length(x) + 2 matrix using  forward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> Array(L̄₁₊(x))\n3×5 Array{Float64,2}:\n 0.0  -1.0   1.0   0.0  0.0\n 0.0   0.0  -1.0   1.0  0.0\n 0.0   0.0   0.0  -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₁₋-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₁₋",
    "category": "method",
    "text": "`L̄₁₋(x)`\n\nReturns a discretized first-order differential operator of length(x) by length(x) + 2 matrix using backward difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> Array(L̄₁₋(x))\n3×5 Array{Float64,2}:\n -1.0   1.0   0.0  0.0  0.0\n  0.0  -1.0   1.0  0.0  0.0\n  0.0   0.0  -1.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L̄₂-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L̄₂",
    "category": "method",
    "text": "`L̄₂(x)`\n\nReturns a discretized second-order differential operator of length(x) by length(x) + 2 matrix using central difference under no boundary condition.\n\nThe first and last columns are applied to the ghost nodes just before x[1] and x[end] respectively.\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> Array(L̄₂(x))\n3×5 Array{Float64,2}:\n 1.0  -2.0   1.0   0.0  0.0\n 0.0   1.0  -2.0   1.0  0.0\n 0.0   0.0   1.0  -2.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊",
    "category": "method",
    "text": "`L₁₊(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`\n\nReturns a discretized first-order differential operator of length(x) by length(x) matrix using forward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper. \n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L₁₊(x, (Reflecting(), Reflecting()))\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0   ⋅\n  0.0  -1.0  1.0\n   ⋅    0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋",
    "category": "method",
    "text": "`L₁₋(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`\n\nReturns a discretized first-order differential operator of length(x) by length(x) matrix using backward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper. \n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L₁₋(x, (Reflecting(), Reflecting()))\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0   ⋅\n -1.0   1.0  0.0\n   ⋅   -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂",
    "category": "method",
    "text": "`L₂(x, bc::Tuple{BoundaryCondition, BoundaryCondition})`\n\nReturns a discretized second-order differential operator of length(x) by length(x) matrix using central difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper. \n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L₂(x, (Reflecting(), Reflecting()))\n3×3 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅\n  1.0  -2.0   1.0\n   ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.x̄-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.x̄",
    "category": "method",
    "text": "`x̄(x)`\n\nReturns an extended grid of length length(x)+2 given grid x.\n\nThe first and last elements of the returned extended grid represent the ghost nodes just before x[1] and x[end] respectively.\n\njulia> x = 1:3\n1:3\n\njulia> x̄(x)\n5-element Array{Int64,1}:\n 0\n 1\n 2\n 3\n 4\n\njulia> x = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> x̄(x)\n5-element Array{Float64,1}:\n 0.5\n 1.0\n 1.5\n 1.7\n 1.9\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}
