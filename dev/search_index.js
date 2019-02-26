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
    "text": "Consider solving for v from the following equation:rho v(x) = f(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v on M-size discretized grids, one can run the following code:# setup \nf(x) = x^2 \nμ = -0.1 # constant negative drift\nσ = 0.1\nM = 100 # size of grid\nx = range(0.0, 1.0, length = M) # grid\n\n# operators with reflecting boundary conditions\nL_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())\nA = μ*L_1_minus + σ^2 / 2 * L_2 \nv_bc = (I * 0.05 - A) \\ f.(x) # solve the value functionNote that the code above uses differential operators with reflecting boundary conditions applied.  One can alternatively use differential operators on interior nodes and stack them with matrices for boundary conditions to compute v:# operators without boundary conditions, adding extra two rows for boundary conditions\nL_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x) # operators on interior nodes\nL = μ*L_1_minus + σ^2 / 2 * L_2 # differential operators on interior nodes\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]]) # matrix for boundary conditions\nA = [([zeros(M) Diagonal(ones(M,M)) zeros(M)] * 0.05 - L); B] # stack them\nv_bar = A \\ [f.(x); 0.0; 0.0] # solve the value function with reflecting barrier bc\nv_interior = v_bar[2:end-1] # extract the interior (is identical with `v_bc`)"
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
    "location": "api/#SimpleDifferentialOperators.diffusionoperators-Tuple{Any,BoundaryCondition,BoundaryCondition}",
    "page": "API",
    "title": "SimpleDifferentialOperators.diffusionoperators",
    "category": "method",
    "text": "diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition)\n\nReturns a tuple of diffusion operators and extended grid (L_1_minus, L_1_plus, L_2, x_bar) with specified boundary conditions.\n\nGiven a grid x of length M, return diffusion operators for negative drift, positive drift, and central differences. BC1 is applied to the lower bound, and BC2 to the upper. x_bar is a (M+2) array that represents the extended grid whose first and last elements represent the ghost nodes just before x[1] and x[end].\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x)\n(L_1_minus =\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [1, 3]  =  0.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [2, 4]  =  0.0\n  [3, 4]  =  1.0, L_1_plus =\n  [1, 2]  =  -1.0\n  [2, 2]  =  0.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -1.0\n  [3, 3]  =  0.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -1.0\n  [3, 5]  =  1.0, L_2 =\n  [1, 1]  =  1.0\n  [1, 2]  =  -2.0\n  [2, 2]  =  1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -2.0\n  [3, 3]  =  1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -2.0\n  [3, 5]  =  1.0, x_bar = [0, 1, 2, 3, 4])\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.diffusionoperators-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.diffusionoperators",
    "category": "method",
    "text": "diffusionoperators(x)\n\nReturns a tuple of diffusion operators and extended grid (L_1_minus, L_1_plus, L_2, x_bar) without applying any boundary conditions.\n\nGiven a grid x of length M, return diffusion operators for negative drift, positive drift, and central differences. No boundary conditions are applied. x_bar is a (M+2) array that represents the extended grid whose first and last elements represent the ghost nodes just before x[1] and x[end].\n\nExamples\n\njulia> x = 1:3\n1:3\n\njulia> L_1_minus, L_1_plus, L_2, x_bar = diffusionoperators(x, Reflecting(), Reflecting())\n(L_1_minus = [0.0 0.0 0.0; -1.0 1.0 0.0; 0.0 -1.0 1.0], L_1_plus = [-1.0 1.0 0.0; 0.0 -1.0 1.0; 0.0 0.0 0.0], L_2 = [-1.0 1.0 0.0; 1.0 -2.0 1.0; 0.0 1.0 -1.0], x_bar = [0, 1, 2, 3, 4])\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}
