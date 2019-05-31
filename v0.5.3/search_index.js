var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#[SimpleDifferentialOperators.jl](https://github.com/QuantEcon/SimpleDifferentialOperators.jl/)-1",
    "page": "Home",
    "title": "﻿SimpleDifferentialOperators.jl",
    "category": "section",
    "text": ""
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "To install, run] add SimpleDifferentialOperatorsNote that this requires Julia 1.1 or later."
},

{
    "location": "#Usage-1",
    "page": "Home",
    "title": "Usage",
    "category": "section",
    "text": "Detailed derivations and more applications can be found here."
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Solving-HJBE-with-constant-drifts-1",
    "page": "Examples",
    "title": "Solving HJBE with constant drifts",
    "category": "section",
    "text": "Consider solving for v from the following equation by the Hamilton-Jacobi-Bellman equation (HJBE):rho v(x) = pi(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x)for some constant rho sigma  0 and mu leq 0. To solve v under the reflecting barrier conditions v(0) = v(1) = 0 on M-size discretized grids, one can run the following code:using LinearAlgebra, SimpleDifferentialOperators\n# setup\nπ(x) = x^2\nμ = -0.1 # constant negative drift\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid (interior points)\n\nx̄ = range(0.0, 1.0, length = (M+2))\nx = interiornodes(x̄) # i.e., x̄[2:end-1]\n\n# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nLₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)\nL_bc = I * ρ - Lₓ\n\n# solve the value function\nv = L_bc \\ π.(x)Note that the interior solution v does not the values of v at the boundary, i.e., v(0) and v(1). To extend the interior solution to the boundary points, one can call extrapolatetoboundary as follows:̄v = extrapolatetoboundary(x̄, v, bc);Here is a complete plot for v:using Plots\nplot(x̄, v̄, lw = 4, label = \"v\")(Image: plot-hjbe-both-reflecting)Note that the code above uses differential operators on the interior nodes with reflecting boundary conditions applied. One can alternatively use operators on extended nodes (extended operators) and stack them with matrices for boundary conditions to compute v:# import SparseArrays package (for identity matrix and spzeros)\nusing SparseArrays\n\n# differential operators on extended nodes\nLₓ = μ*L₁₋(x̄) + σ^2 / 2 * L₂(x̄)\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[-1; 1; zeros(M)] [zeros(M); -1; 1]])\nb = [0.0; 0.0]\n\n# form bellman equation on extension\nL = [spzeros(M) ρ*I spzeros(M)] - Lₓ\n\n# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L; B] \\ [π.(x); b]\n\n# extract the interior (is identical with `v` above)\nv =  v̄[2:end-1]"
},

{
    "location": "examples/#Solving-HJBE-with-absorbing-barrier-conditions-1",
    "page": "Examples",
    "title": "Solving HJBE with absorbing barrier conditions",
    "category": "section",
    "text": "Instead of having the reflecting barrier conditions on both lower bound and upper bound v(0) = v(1) = 0 as above, one can impose an absorbing barrier condition as well. To solve v under the reflecting barrier conditions v(0) = S (absorbing barrier on lower bound) for some S and v(1) = 0 (reflecting barrier on upper bound), one can construct B and b for the boundary conditions as follows.First, consider the case where S neq 0, which gives a nonhomogenous boundary condition:# define S\nS = 3.0\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])\nb = [S; 0.0];We can then apply one Gaussian elimination step to remove a non-zero element of the first column in L, which is mu Delta^-1 - (sigma^22) Delta^-2. This can be done by substracting the first row of the stacked system L B by the first row of the system B = b by mu Delta^-1 - (sigma^22) Delta^-2. This returns the following identical system:beginbmatrix\nL2M+1 \nB2\nendbmatrix\n=\nbeginbmatrix\nπ^* \nb2\nendbmatrixwhereπ^* =\nbeginbmatrix\nπ(x_1) - S(mu Delta^-1 - (sigma^22) Delta^-2)\n \nvdots\n\nπ(x_M)\nendbmatrixNow solve v:# stack the systems of bellman and boundary conditions, and solve\nv̄ =  [L; B] \\ [π.(x); b]Here is a plot for v:plot(x̄, v̄, lw = 4, label = \"v\")(Image: plot-hjbe-lb-absorbing-ub-reflecting)Note that this can be alternatively done by constructing the corresponding differential operators on the interior with Absorbing() boundary condition when S = 0:# discretize L = ρ - μ D_x - σ^2 / 2 D_xx\n# subject to reflecting barriers at 0 and 1\nbc = (Absorbing(), Reflecting())\nLₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄ , bc)\nL_bc = I * ρ - Lₓ\n\n# solve the value function\nv = L_bc \\ π.(x)In fact, on the interior, they return identical solutions:# define S\nS = 0.0\n\n# boundary conditions (i.e. B v̄ = b)\nB = transpose([[1; 0; zeros(M)] [zeros(M); -1; 1]])\nb = [S; 0.0];\n\n# stack the systems of bellman and boundary conditions, and solve\nv̄ = [L; B] \\ [π.(x); b]\n\n# confirm that v returns the identical solution as the one from the stacked system\nusing Test\n@test v ≈ v̄[2:end-1]"
},

{
    "location": "examples/#Solving-HJBE-with-jump-diffusion-1",
    "page": "Examples",
    "title": "Solving HJBE with jump diffusion",
    "category": "section",
    "text": "Consider the jump process added to the HJBE with some intensity lambda geq 0: rho v(x) = pi(x) + mu partial_x v(x) + fracsigma^22 partial_xx v(x) + lambda left v(x + J(x) ) - v(x) rightwhere the jump process is defined by the jump magnitude defined by J(x_i). In SimpleDifferentialOperators.jl, the jump process can be defined as follows:# uniform jump\njumpf(x_i) = -0.01\njumpprocess = JumpProcess(x̄, jumpf)Note that, the corresponding indices for destinations will be determined by the nearest neighbor as the domain has to be discretized accordingly. Alternatively, if the jump magnitude is uniform across all cohorts, one can forward the uniform jump magnitude as follows:# use the fact that the jump magnitude is uniform across all nodes\njumpprocess = JumpProcess(x̄, -0.01)One can define a jump process manually by providing jump magnitudes in indices as well. If a jump process is defined by the indices on a discretized domain, incurring jumps from v(x_i) to v(x_i-1) for all i in 2 leq i leq M, one can construct a jump process as follows:# length of nodes on the interior\nM = length(interiornodes(x̄))\n# vector of jumps; ith element represents the jump from ith node in the interior\njumps = -ones(M)\n# define jump process \njumpprocess =  JumpProcess(x̄, jumps)Alternatively, one can define an identical jump process with ease if the jump maginitude in index is uniform across all nodes:# use the fact that the jump size is uniform across all nodes\njumpprocess = JumpProcess(x̄, -1)Then one can define the corresponding discretized operator L_n and solve value functions as follows:# define jump intensity\nλ = 0.6\n\n# construct discretized operators on interior nodes with the bc\nLₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄ , bc) + λ * Lₙbc(x̄, bc, jumpprocess) \nL_bc = I * ρ - Lₓ\n\n# solve the value function\nv = L_bc \\ π.(x)"
},

{
    "location": "examples/#Solving-HJBE-with-jump-diffusions-and-Markov-chains-1",
    "page": "Examples",
    "title": "Solving HJBE with jump diffusions and Markov chains",
    "category": "section",
    "text": "Suppose we are asked to solve HJBE with two states (N=2) where for each ith state with the corresponding differential operator L_i under different payoff functions pi_i and drifts mu_i, there is a transition intensity of q_ij to have state j assigned.# setup\n# payoff functions\nπ_1(x) = x^2\nπ_2(x) = (x-0.01)^2\n\n# constant negative drifts\nμ_1 = -0.1\nμ_2 = -0.15\nλ = 0.6\nσ = 0.1\nρ = 0.05\nM = 100 # size of grid (interior points)\n\nx̄ = range(0.0, 1.0, length = (M+2))\nx = interiornodes(x̄) # i.e., x̄[2:end-1]Let the HJBE in the first state have a jump process J associated while the one for the second state does not. Then we have the following system of differential equations; note that we have q_ii = -q_ij for i neq j:beginalign\nrho v_1 (x) = pi_1(x) + mu_1 partial_x v_1(x) + fracsigma^22 partial_xx v_1(x) + lambda left v_1(x + J(x) ) - v_1(x) right + q_12  v_2(x) - v_1(x)  \nrho v_2 (x) = pi_2(x) + mu_2 partial_x v_2(x) + fracsigma^22 partial_xx v_2(x) + q_21  v_1(x) - v_2(x) \nendalignFirst, construct L_1 and L_2 ignoring the Markov chain for transition between the two states; assume that both states have reflecting boundary conditions applied.# construct the jump process for the operator in state 1\njumpprocess1 = JumpProcess(x̄, -0.01)\n\n# construct the differential operators for both states\n# subject to reflecting barriers at 0 and 1\nbc = (Reflecting(), Reflecting())\nL_1ₓ = μ_1*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc) + λ * Lₙbc(x̄, bc, jumpprocess1)\nL_2ₓ = μ_2*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)\nL_1_bc = I * ρ - L_1ₓ\nL_2_bc = I * ρ - L_2ₓThen construct an intensity matrix Q, whose (ij)th element represents q_ij:# define intensity matrix for transition\nQ = [-0.01 0.01; 0.02 -0.02]Using the discretized operators L_1_bc and L_2_bc on interior nodes x with boundary conditions bc applied, one can construct the joint operator L_bc with the intensity matrix Q as follows:# define the corresponding joint operator\nL_bc = jointoperator_bc((L_1_bc, L_2_bc), Q)Construct a vector of payoff functions pi_1 and pi_2 stacked together and solve the system using the joint operator constructed above:# solve the system\nv = L_bc \\ [π_1.(x); π_2.(x)]Note that the first M elements represent the discretized solution for v_1 and the last M elements represent the one for v_2:# extract the solution for each state\nv_1 = v[1:M]\nv_2 = v[(M+1):end]\n\n# plot v_1 and v_2\nplot(x, [v_1, v_2], lw = 4, label = [\"v_1\", \"v_2\"])(Image: plot-hjbe-two-states)"
},

{
    "location": "examples/#Solving-HJBE-with-state-dependent-drifts-1",
    "page": "Examples",
    "title": "Solving HJBE with state-dependent drifts",
    "category": "section",
    "text": "One can also deploy upwind schemes when drift variable is not constant. Consider solving for v from the following Bellman equation:rho v(x) = π(x) + mu(x) partial_x v(x) + fracsigma^22 partial_xx v(x)associated with the diffusion processdx = mu(x) dt + sigma dWfor some constant rho sigma  0 and mu(x) = -x. Note that mu(x) depends on states. The following code will solve v using upwind schemes, with the reflecting barrier conditions v(0) = v(1) = 0 applied:# setup\nπ(x) = x^2\nμ(x) = -x # drift depends on state\nσ = 1.0\nρ = 0.05\nM = 100 # size of grid\n\nx̄ = range(-1., 1., length = M + 2)\nx = interiornodes(x̄) # i.e., x̄[2:end-1]\n\nbc = (Reflecting(), Reflecting())\n\n# Define first order differential operator using upwind scheme\nL₁ = Diagonal(min.(μ.(x), 0.0)) * L₁₋bc(x̄, bc) + Diagonal(max.(μ.(x), 0.0)) * L₁₊bc(x̄, bc)\n\n# Define linear operator using upwind schemes\nLₓ = L₁ - σ^2 / 2 * L₂bc(x̄, bc)\nL_bc_state_dependent = I * ρ - Lₓ\n\n# solve the value function\nv = L_bc_state_dependent \\ π.(x)"
},

{
    "location": "examples/#Finding-stationary-distribution-from-the-Kolmogorov-forward-equation-(KFE)-1",
    "page": "Examples",
    "title": "Finding stationary distribution from the Kolmogorov forward equation (KFE)",
    "category": "section",
    "text": "The KFE equation ispartial_t f(xt) = -mu partial_x f(xt) + fracsigma^22 partial_xx f(xt)for x in (x_min x_max) with the following corresponding reflecting barrier conditions:beginalign\n-mu f(x_min t) +fracsigma^22 partial_x f(x_min t) = 0 \n-mu f(x_max t) +fracsigma^22 partial_x f(x_max t) = 0\nendaligni.e.,beginalign\n-frac2musigma^2 f(x_min t) +partial_x f(x_min t) = 0 \n-frac2musigma^2 f(x_max t) +partial_x f(x_max t) = 0\nendalignwhich gives mixed boundary conditions with overlinexi = underlinexi = -frac2musigma^2.One can compute the stationary distribution of the state x above from the corresponding KFE by taking partial_t f(xt) = 0, i.e., solving f from the L^* f(x) = 0 whereL^* = - mu partial_x + fracsigma^22 partial_xxThe following code constructs L^*:# parameter setup\nμ = -0.1 # constant negative drift\nσ = 0.1\nM = 100 # size of grid (interior points)\nx_min = 0.0\nx_max = 1.0\nx̄ = range(x_min, x_max, length = (M+2))\n\n# ξ values for mixed boundary conditions\nξ_lb = ξ_ub = -2μ/σ^2\n\n# define the corresponding mixed boundary conditions\n# note that the direction on the lower bound is backward (default is forward)\n# as the drift μ is negative.\nbc = (Mixed(ξ = ξ_lb, direction = :backward), Mixed(ξ = ξ_ub))\n\n# use SimpleDifferentialOperators.jl to construct the operator on the interior\nL_KFE = Array(-μ*L₁₊bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc))One can find the stationary distribution f by solving the following discretized system of equations:L^* f = 0such that the sum of f is one. This can be found by finding a non-trivial eigenvector f_ss for L_KFE associated with the eigenvalue of zero:using Arpack # library for extracting eigenvalues and eigenvectors\n\n# extract eigenvalues and eigenvectors, smallest eigenval in magintute first\nλ, ϕ = eigs(L_KFE, which = :SM); \n# extract the very first eigenvector (associated with the smallest eigenvalue)\nf_ss = real.(ϕ[:,1]);\n# normalize it\nf_ss = f_ss / sum(f_ss)Using L from the state-dependent drift example above, this results in the following stationary distribution:plot(x, f_ss, lw = 4, label = \"f_ss\")(Image: plot-stationary-dist)Note that the operator for the KFE in the original equation is the adjoint of the operator for infinitesimal generator used in the HJBE, L, and the correct discretization scheme for L^* is, analogously, done by taking the transpose of the discretized operator for HJBE, L (See Gabaix et al., 2016 and Achdou et al., 2017), which has been constructed as Lₓ is the HJBE example above. In fact, the discretized L^* and L^T are identical:# discretize L = μ D_x + σ^2 / 2 D_xx\n# for infinitesimal generators used in the HJBE\n# subject to reflecting barrier conditions\nbc = (Reflecting(), Reflecting())\nLₓ = μ*L₁₋bc(x̄, bc) + σ^2 / 2 * L₂bc(x̄, bc)\n\n@test transpose(Lₓ) == L_KFE"
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
    "text": "We prepared a gallery of notebooks to showcase the applicability of the package to different economic problems.Kolmogorov Forward Equations: HTML and Jupyter notebook\nOptimal Stopping Problems: HTML and Jupyter notebookOptimal Stopping Problems with Jump Diffusion: HTML and Jupyter notebook\nComputational Appendix for Optimal Stopping Problems: HTML and Jupyter notebook"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#SimpleDifferentialOperators.DifferentialOperator-Tuple{Any,Tuple{BoundaryCondition,BoundaryCondition},BackwardFirstDifference}",
    "page": "API",
    "title": "SimpleDifferentialOperators.DifferentialOperator",
    "category": "method",
    "text": "DifferentialOperator(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition}, method::DiscretizationMethod)\n\nReturns a discretized differential operator of  length(interiornodes(x̄)) by length(interiornodes(x̄)) matrix under mixed boundary conditions from bc using a discretization method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), BackwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0    ⋅    ⋅\n -1.0   1.0   0.0   ⋅\n   ⋅   -1.0   1.0  0.0\n   ⋅     ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), ForwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅    ⋅\n  0.0  -1.0   1.0   ⋅\n   ⋅    0.0  -1.0  1.0\n   ⋅     ⋅    0.0  0.0\n\njulia> DifferentialOperator(x̄, (Reflecting(), Reflecting()), CentralSecondDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  1.0  -2.0   1.0    ⋅\n   ⋅    1.0  -2.0   1.0\n   ⋅     ⋅    1.0  -1.0\n\njulia> x̄ = 0:5\n0:5\n\njulia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), BackwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n Inf     0.0    ⋅    ⋅\n  -1.0   1.0   0.0   ⋅\n    ⋅   -1.0   1.0  0.0\n    ⋅     ⋅   -1.0  1.0\n\njulia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), ForwardFirstDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  0.0  -1.0   1.0    ⋅\n   ⋅    0.0  -1.0   1.0\n   ⋅     ⋅    0.0  -0.5\n\njulia> DifferentialOperator(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)), CentralSecondDifference())\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -Inf     1.0    ⋅     ⋅\n    1.0  -2.0   1.0    ⋅\n     ⋅    1.0  -2.0   1.0\n     ⋅     ⋅    1.0  -1.5\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.ExtensionDifferentialOperator-Tuple{AbstractRange,BackwardFirstDifference}",
    "page": "API",
    "title": "SimpleDifferentialOperators.ExtensionDifferentialOperator",
    "category": "method",
    "text": "ExtensionDifferentialOperator(x̄, method::DiscretizationMethod)\n\nReturns a discretized differential operator of  length(interiornodes(x̄)) by length(x̄) matrix under no boundary condition using a discretization method specified by method.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> ExtensionDifferentialOperator(x̄, BackwardFirstDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 1]  =  -1.0\n  [1, 2]  =  1.0\n  [2, 2]  =  -1.0\n  [2, 3]  =  1.0\n  [3, 3]  =  -1.0\n  [3, 4]  =  1.0\n  [4, 4]  =  -1.0\n  [4, 5]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, ForwardFirstDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 8 stored entries:\n  [1, 2]  =  -1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -1.0\n  [3, 5]  =  1.0\n  [4, 5]  =  -1.0\n  [4, 6]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, CentralSecondDifference())\n4×6 SparseArrays.SparseMatrixCSC{Float64,Int64} with 12 stored entries:\n  [1, 1]  =  1.0\n  [1, 2]  =  -2.0\n  [2, 2]  =  1.0\n  [1, 3]  =  1.0\n  [2, 3]  =  -2.0\n  [3, 3]  =  1.0\n  [2, 4]  =  1.0\n  [3, 4]  =  -2.0\n  [4, 4]  =  1.0\n  [3, 5]  =  1.0\n  [4, 5]  =  -2.0\n  [4, 6]  =  1.0\n\njulia> ExtensionDifferentialOperator(x̄, JumpProcess(x̄, -1.0))\n4×6 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:\n 0.0  0.0    ⋅     ⋅     ⋅    ⋅\n  ⋅   1.0  -1.0    ⋅     ⋅    ⋅\n  ⋅    ⋅    1.0  -1.0    ⋅    ⋅\n  ⋅    ⋅     ⋅    1.0  -1.0   ⋅\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₀affine-Tuple{Any,Tuple{BoundaryCondition,BoundaryCondition}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₀affine",
    "category": "method",
    "text": "L₁₋affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`.\n\nReturns a `length(interiornodes(x̄))`-length vector such that solving\n\nL₁₋ * v(x̄) = b on interior nodes x under the boundary condition given by bc is identical with L₁₋bc * v(x) = b + L₁₋affine(x̄, bc) The first element of bc is  applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₀affine(x̄, (Reflecting(), Reflecting()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> L₀affine(x̄, (NonhomogeneousAbsorbing(1.0), Reflecting()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> L₀affine(x̄, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊",
    "category": "method",
    "text": "L₁₊(x̄)\n\nReturns a discretized first-order differential operator of length(interiornodes(x̄)) by length(x̄) matrix using forward difference under no boundary condition.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> Array(L₁₊(x̄))\n4×6 Array{Float64,2}:\n 0.0  -1.0   1.0   0.0   0.0  0.0\n 0.0   0.0  -1.0   1.0   0.0  0.0\n 0.0   0.0   0.0  -1.0   1.0  0.0\n 0.0   0.0   0.0   0.0  -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊affine-Tuple{Any,Tuple{BoundaryCondition,BoundaryCondition}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊affine",
    "category": "method",
    "text": "L₁₊affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`.\n\nReturns a `length(interiornodes(x̄))`-length vector such that solving\n\nL₁₋ * v(x̄) = b on interior nodes x under the boundary condition given by bc is identical with L₁₋bc * v(x) = b + L₁₊affine(x̄, bc) The first element of bc is  applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₊affine(x̄, (Reflecting(), Reflecting()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> L₁₊affine(x̄, (NonhomogeneousAbsorbing(1.0), Reflecting()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> L₁₊affine(x̄, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))\n4-element Array{Float64,1}:\n  0.0\n  0.0\n  0.0\n -2.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₊bc-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₊bc",
    "category": "method",
    "text": "L₁₊bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of  length(interiornodes(x̄)) by length(interiornodes(x̄)) matrix using forward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₊bc(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅    ⋅\n  0.0  -1.0   1.0   ⋅\n   ⋅    0.0  -1.0  1.0\n   ⋅     ⋅    0.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋",
    "category": "method",
    "text": "L₁₋(x̄)\n\nReturns a discretized first-order differential operator of  length(interiornodes(x̄)) by length(x̄) matrix using backward difference under no boundary condition.\n\nExamples\n\njulia> x̄ = 1:3\n1:3\n\njulia> Array(L₁₋(x̄))\n1×3 Array{Float64,2}:\n -1.0  1.0  0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋affine-Tuple{Any,Tuple{BoundaryCondition,BoundaryCondition}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋affine",
    "category": "method",
    "text": "L₁₋affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`.\n\nReturns a `length(interiornodes(x̄))`-length vector such that solving\n\nL₁₋ * v(x̄) = b on interior nodes x under the boundary condition given by bc is identical with L₁₋bc * v(x) = b + L₁₋affine(x̄, bc). The first element of bc is  applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₋affine(x̄, (Reflecting(), Reflecting()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> L₁₋affine(x̄, (NonhomogeneousAbsorbing(1.0), Reflecting()))\n4-element Array{Float64,1}:\n 1.0\n 0.0\n 0.0\n 0.0\n\njulia> L₁₋affine(x̄, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))\n4-element Array{Float64,1}:\n 1.0\n 0.0\n 0.0\n 0.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₁₋bc-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₁₋bc",
    "category": "method",
    "text": "L₁₋bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized first-order differential operator of  length(interiornodes(x̄)) by length(interiornodes(x̄)) matrix using backward difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₁₋bc(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0    ⋅    ⋅\n -1.0   1.0   0.0   ⋅\n   ⋅   -1.0   1.0  0.0\n   ⋅     ⋅   -1.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂",
    "category": "method",
    "text": "L₂(x̄)\n\nReturns a discretized second-order differential operator of length(interiornodes(x̄)) by length(x̄) matrix using central difference under no boundary condition.\n\nExamples\n\njulia> x̄ = 0:5\n0:5 \n\njulia> Array(L₂(x̄))\n4×6 Array{Float64,2}:\n 1.0  -2.0   1.0   0.0   0.0  0.0\n 0.0   1.0  -2.0   1.0   0.0  0.0\n 0.0   0.0   1.0  -2.0   1.0  0.0\n 0.0   0.0   0.0   1.0  -2.0  1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂affine-Tuple{Any,Tuple{BoundaryCondition,BoundaryCondition}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂affine",
    "category": "method",
    "text": "L₂affine(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a vector for affine operator given an extended grid `x̄` and boundary conditions `bc`.\n\nReturns a `length(interiornodes(x̄))`-length vector such that solving\n\nL₂ * v(x̄) = b on interior nodes x under the boundary condition given by bc is identical with L₂bc * v(x) = b + L₂affine(x̄, bc) The first element of bc is  applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₂affine(x̄, (Reflecting(), Reflecting()))\n4-element Array{Float64,1}:\n 0.0\n 0.0\n 0.0\n 0.0\n\njulia> L₂affine(x̄, (NonhomogeneousAbsorbing(1.0), Reflecting()))\n4-element Array{Float64,1}:\n -1.0\n  0.0\n  0.0\n  0.0\n\njulia> L₂affine(x̄, (NonhomogeneousAbsorbing(1.0), NonhomogeneousAbsorbing(2.0)))\n4-element Array{Float64,1}:\n -1.0\n  0.0\n  0.0\n -2.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.L₂bc-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.L₂bc",
    "category": "method",
    "text": "L₂bc(x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a discretized second-order differential operator of  length(interiornodes(x̄)) by length(interiornodes(x̄)) matrix using central difference under boundary conditions specified by bc.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L₂bc(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  1.0  -2.0   1.0    ⋅\n   ⋅    1.0  -2.0   1.0\n   ⋅     ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.Lₙ-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.Lₙ",
    "category": "method",
    "text": "Lₙ(x̄, method)\n\nReturns a discretized jump process operator of length(interiornodes(x̄)) by length(x̄) matrix specified by method\n\nExamples\n\njulia> x̄ = 0:5\n0:5 \n\njulia> Lₙ(x̄, JumpProcess(x̄, -1.0))\n4×6 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:\n 0.0  0.0    ⋅     ⋅     ⋅    ⋅\n  ⋅   1.0  -1.0    ⋅     ⋅    ⋅\n  ⋅    ⋅    1.0  -1.0    ⋅    ⋅\n  ⋅    ⋅     ⋅    1.0  -1.0   ⋅\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.Lₙbc-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.Lₙbc",
    "category": "method",
    "text": "Lₙbc(x̄, method, (Absorbing(), Absorbing()))\n\nReturns a discretized jump process operator of length(interiornodes(x̄)) by length(interiornodes(x̄)) matrix specified by method under boundary conditions specified by bc.\n\nExamples\n\njulia> x̄ = 0:5\n0:5 \n\njulia> Lₙbc(x̄, (Absorbing(), Absorbing()), JumpProcess(x̄, -1.0))\n4×4 BandedMatrices.BandedMatrix{Float64,Array{Float64,2},Base.OneTo{Int64}}:\n 0.0   0.0    ⋅     ⋅\n 1.0  -1.0   0.0    ⋅\n  ⋅    1.0  -1.0   0.0\n  ⋅     ⋅    1.0  -1.0\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.extrapolatetoboundary-Tuple{Any,Any,Tuple{BoundaryCondition,BoundaryCondition}}",
    "page": "API",
    "title": "SimpleDifferentialOperators.extrapolatetoboundary",
    "category": "method",
    "text": "extrapolatetoboundary(v, x̄, bc::Tuple{BoundaryCondition, BoundaryCondition})\n\nReturns a length(x̄)-vector whose 2:(length(x̄)-1) elements are v, the first and last element are extrapolated v on the boundaries of x̄ according to boundary conditions bc given.\n\nThe first element of bc is applied to the lower bound, and second element of bc to the upper.\n\nExamples\n\njulia> x̄ = -2:2\n-2:2\n\njulia> x = interiornodes(x̄)\n-1:1\n\njulia> v = (x -> x^2).(x)\n3-element Array{Int64,1}:\n 1\n 0\n 1\n\njulia> extrapolatetoboundary(v, x̄, (Absorbing(), Absorbing()))\n5-element Array{Int64,1}:\n 0\n 1\n 0\n 1\n 0\n\njulia> extrapolatetoboundary(v, x̄, (Absorbing(), Reflecting()))\n5-element Array{Int64,1}:\n 0\n 1\n 0\n 1\n 1\n\njulia> extrapolatetoboundary(v, x̄, (Mixed(ξ = 3.0), Reflecting()))\n5-element Array{Float64,1}:\n -0.5\n  1.0\n  0.0\n  1.0\n  1.0\n  \n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.interiornodes-Tuple{Any,Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.interiornodes",
    "category": "method",
    "text": "interiornodes(x̄, bc)\n\nReturns an interior grid corresponding to the boundary condition bc given extended grid x̄.\n\njulia> x̄ = 0:5\n0:5\n\njulia> interiornodes(x̄, (Reflecting(), Reflecting()))\n1:4\n\njulia> x̄ = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> interiornodes(x̄, (Mixed(ξ = 1.0), Mixed(ξ = 1.0)))\n1-element Array{Float64,1}:\n 1.5\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.interiornodes-Tuple{Any}",
    "page": "API",
    "title": "SimpleDifferentialOperators.interiornodes",
    "category": "method",
    "text": "interiornodes(x̄)\n\nReturns an interior grid of length length(x̄)-2 given extended grid x̄.\n\njulia> x̄ = 0:5\n0:5\n\njulia> interiornodes(x̄)\n1:4\n\njulia> x̄ = [1.0; 1.5; 1.7]\n3-element Array{Float64,1}:\n 1.0\n 1.5\n 1.7\n\njulia> interiornodes(x̄)\n1-element Array{Float64,1}:\n 1.5\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.jointoperator_bc-Tuple{Any,Array}",
    "page": "API",
    "title": "SimpleDifferentialOperators.jointoperator_bc",
    "category": "method",
    "text": "jointoperator_bc(operators, Q::Array)\n\nReturns a discretized operator that solves systems of differential equations defined by operators with transitions by Q where operators is an N-length collection of  discretized operators with boundary conditions applied and Q is N by N intensity matrix.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> L1bc = L₁₋bc(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n  0.0   0.0    ⋅    ⋅\n -1.0   1.0   0.0   ⋅\n   ⋅   -1.0   1.0  0.0\n   ⋅     ⋅   -1.0  1.0\n\njulia> L2bc = L₂bc(x̄, (Reflecting(), Reflecting()))\n4×4 LinearAlgebra.Tridiagonal{Float64,Array{Float64,1}}:\n -1.0   1.0    ⋅     ⋅\n  1.0  -2.0   1.0    ⋅\n   ⋅    1.0  -2.0   1.0\n   ⋅     ⋅    1.0  -1.0\n\njulia> Q = [-0.5 0.5; 0.3 -0.3]\n2×2 Array{Float64,2}:\n -0.5   0.5\n  0.3  -0.3\n\njulia> jointoperator_bc((L1bc, L2bc), Q)\n2×2-blocked 8×8 BlockBandedMatrices.BandedBlockBandedMatrix{Float64,BlockArrays.PseudoBlockArray{Float64,2,Array{Float64,2},BlockArrays.BlockSizes{2,Tuple{Array{Int64,1},Array{Int64,1}}}}}:\n -0.5   0.0    ⋅    ⋅   │   0.5   0.0    ⋅     ⋅\n -1.0   0.5   0.0   ⋅   │   0.0   0.5   0.0    ⋅\n   ⋅   -1.0   0.5  0.0  │    ⋅    0.0   0.5   0.0\n   ⋅     ⋅   -1.0  0.5  │    ⋅     ⋅    0.0   0.5\n ───────────────────────┼────────────────────────\n  0.3   0.0    ⋅    ⋅   │  -1.3   1.0    ⋅     ⋅\n  0.0   0.3   0.0   ⋅   │   1.0  -2.3   1.0    ⋅\n   ⋅    0.0   0.3  0.0  │    ⋅    1.0  -2.3   1.0\n   ⋅     ⋅    0.0  0.3  │    ⋅     ⋅    1.0  -1.3\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.JumpProcess",
    "page": "API",
    "title": "SimpleDifferentialOperators.JumpProcess",
    "category": "type",
    "text": "JumpProcess(x̄, jumps::AbstractArray, truncate = (:interior, :interior))\n\nReturns a DiscretizationMethod object that can be used to construct  a discretized operator jump process. \n\njumps is a (length(x̄)-2)-vector whose ith element is an integer that represents  a jump size in index and direction (by sign) from ith element of interiornodes(x̄) and the first and second elements of truncate represent truncation location for  the lower bound and upper bound when the jump is out of the truncated boundary  (:interior for the first/last element of interiornodes(̄x) and :exterior for  the first/last element of x̄). The default parameter is (:interior, :interior).\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> JumpProcess(x̄, [-1; -1; -1; -1], (:interior, :interior))\nJumpProcess{Array{Int64,1}}([0, -1, -1, -1])\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.JumpProcess",
    "page": "API",
    "title": "SimpleDifferentialOperators.JumpProcess",
    "category": "type",
    "text": "JumpProcess(x̄, jumpf::Function, truncate = (:interior, :interior))\n\nReturns a DiscretizationMethod object that can be used to construct  a discretized operator jump process. \n\njumpf is a function that takes an element of interiornodes(x̄) and returns  a scalar that represents the corresponding nominal jump size and direction (by sign) and the first and second elements of truncate represent truncation location for  the lower bound and upper bound when the jump is out of the truncated boundary  (:interior for the first/last element of interiornodes(̄x) and :exterior for  the first/last element of x̄). The default parameter is (:interior, :interior). The code uses nearest-neighbour rule to determine the indices of destinations according to `jumpf.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> jumpf(x) = -1.4\njumpf (generic function with 1 method)\n\njulia> JumpProcess(x̄, jumpf)\nJumpProcess{Array{Int64,1}}([0, -1, -1, -1])\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.JumpProcess",
    "page": "API",
    "title": "SimpleDifferentialOperators.JumpProcess",
    "category": "type",
    "text": "JumpProcess(x̄, uniformjumpsize::Real, truncate = (:interior, :interior))\n\nReturns a DiscretizationMethod object that can be used to construct  a discretized operator jump process. \n\nuniform_jump_size is a scalar in Real that represents a nominal jump size  and direction (by sign) from all elements of interiornodes(x̄) and the first and second elements of truncate represent truncation location for  the lower bound and upper bound when the jump is out of the truncated boundary  (:interior for the first/last element of interiornodes(̄x) and :exterior for  the first/last element of x̄). The default parameter is (:interior, :interior). The code uses nearest-neighbour rule to determine the indices of destinations according to `jumpf.\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> JumpProcess(x̄, -1.4)\nJumpProcess{Array{Int64,1}}([0, -1, -1, -1])\n\n\n\n\n\n"
},

{
    "location": "api/#SimpleDifferentialOperators.JumpProcess",
    "page": "API",
    "title": "SimpleDifferentialOperators.JumpProcess",
    "category": "type",
    "text": "JumpProcess(x̄, uniform_jump::AbstractArray, truncate = (:interior, :interior))\n\nReturns a DiscretizationMethod object that can be used to construct  a discretized operator jump process. \n\nuniform_jump is a scalar Int64 that represents a jump size in index  and direction (by sign) from all elements of interiornodes(x̄) and the first and second elements of truncate represent truncation location for  the lower bound and upper bound when the jump is out of the truncated boundary  (:interior for the first/last element of interiornodes(̄x) and :exterior for  the first/last element of x̄). The default parameter is (:interior, :interior).\n\nExamples\n\njulia> x̄ = 0:5\n0:5\n\njulia> JumpProcess(x̄, -1)\nJumpProcess{Array{Int64,1}}([0, -1, -1, -1])\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}