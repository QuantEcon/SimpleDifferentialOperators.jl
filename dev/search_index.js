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
    "text": "To install, run] dev https://github.com/QuantEcon/SimpleDifferentialOperators.jl"
},

{
    "location": "#Usage-1",
    "page": "﻿Installation",
    "title": "Usage",
    "category": "section",
    "text": "robin_diffusionoperators(x, ξ) takes a vector or range of grids x and a coefficient ξ to be used for boundary conditions\nbeginaligned\n    xi v(0) + nabla_z v(0 ) = 0 \n    xi v(barz) + nabla_z v(barz) = 0\nendaligned\n\n \nand returns a tuple of discretized differential operators in matrix (L_1_minus L_1_plus L_2) that represent the corresponding forward-difference first-order derivative operator backward-difference first-order derivative operator and the centralized second-order derivative operator respectively \n\nFor instance to compute the first-order derivative of scalar-valued function f on given grids x with the boundary condition coefficient ξ one can run the following lines\n\njulia\nL_1_plus = robin_diffusionoperator(x ξ)  compute discretized differential operator\nf_prime_fd = L_1_plus * f(x)  first order derivative using forward difference\n\n\nWhen ξ is zero reflecting_diffusionoperators(x) can be used instead\n\n\nDerivation\n=====\nMore details on derivation can be found here(httpsgithubcomubceconcomputing_and_datascienceblobmastercontinuous_time_methodsnotesdifferential-operator-on-irregular-gridstex)\n\nSetup\n-------------\n\n-   Define an irregular grid z_i_i=1^P with z_1 = 0 and\n    z_P = barz is a large number Denote the grid with the\n    variable name ie z equiv z_i_i=1^P\n\n-   Denote the distance between the grid points as the *backward*\n    difference\n\n    \n    beginaligned\n        Delta_i- equiv z_i - z_i-1 textfor  i = 2ldots P\n        Delta_i+ equiv z_i+1 - z_i textfor  i = 1ldots P-1\n    endaligned\n    \n\n-   Assume Delta_1 - = Delta_1 + and\n    Delta_P + = Delta_P - due to ghost points z_0 and\n    z_P+1 on both boundaries (ie the distance to the ghost nodes\n    are the same as the distance to the closest nodes) Then define the\n    vector of backward and forward first differences as \n\n    \n        Delta_- equiv beginbmatrix z_2 - z_1 \n        textdiff(z)\n        endbmatrix\n    \n\n    \n        Delta_+ equiv beginbmatrix textdiff(z)\n        z_P - z_P-1\n        endbmatrix\n    \n    \n\n-   Reflecting barrier conditions \n\n    \n    beginaligned\n        xi v(0) + nabla_z v(0 ) = 0  (1) \n        xi v(barz) + nabla_z v(barz) = 0  (2)\n    endaligned\n    \n    \nLet L_1^- be the discretized backward first differences and L_2\nbe the discretized central differences subject to the Neumann boundary\nconditions such that L_1^- v(t) and L_2 v(t) represent the first\nand second derivatives of v(z) respectively at t For second\nderivatives we use the following numerical scheme from Achdou et al\n(2017)(httpwwwprincetonedumollHACTpdf)\n\nv(z_i) approx dfrac Delta_i- v( z_i + Delta_i+) - (Delta_i+ + Delta_i-) v( z_i ) + Delta_i+ v( z_i - Delta_i-)frac12(Delta_i+ + Delta_i-) Delta_i+ Delta_i-  textfor  i = 1 ldots P"
},

{
    "location": "#Regular-grids-1",
    "page": "﻿Installation",
    "title": "Regular grids",
    "category": "section",
    "text": "Suppose that the grids are regular, i.e., elements of textdiff(z) are all identical with Delta for some Delta  0.Using the backward first-order difference, the first boundary condition (1) can be alternatively represented as $\\begin{aligned} \\dfrac{v(0) - v(-\\Delta)}{\\Delta} &= - \\xi v(0) \\end{aligned} $Similarly, using discretized central differences of second orders, the first boundary condition (1) yields$\\begin{aligned} \\dfrac{v (\\Delta) - 2 v(0) + v(-\\Delta)}{\\Delta^2} &=   \\dfrac{v(\\Delta) - v(0)}{\\Delta^2} - \\dfrac{1}{\\Delta}\\dfrac{v (0) - v(-\\Delta) }{\\Delta}  \\\n&= \\dfrac{v(\\Delta) - v(0)}{\\Delta^2} + \\dfrac{1}{\\Delta} \\xi v(0)  \\  &= \\dfrac{1}{\\Delta^2}  (- 1 + \\Delta \\xi) v(0)  + \\dfrac{1}{\\Delta^2}  v(\\Delta)   \\end{aligned} $Similarly, for the secondary boundary condition (2), we have $\\begin{aligned} \\dfrac{v (\\bar{z} + \\Delta) - 2 v(\\bar{z} ) + v(\\bar{z} -\\Delta)}{\\Delta^2} &=   \\dfrac{v(\\bar{z} - \\Delta) - v(\\bar{z})}{\\Delta^2} + \\dfrac{1}{\\Delta}\\dfrac{ v(\\bar{z}+\\Delta) - v (\\bar{z}) }{\\Delta}  \\\n&= \\dfrac{v(\\bar{z} - \\Delta) - v(\\bar{z})}{\\Delta^2}  - \\dfrac{1}{\\Delta} \\xi v(\\bar{z})  \\  &= \\dfrac{1}{\\Delta^2}  (- 1 - \\Delta \\xi) v(\\bar{z})  + \\dfrac{1}{\\Delta^2}  v(\\bar{z} - \\Delta)   \\end{aligned} $Thus, the corresponding L_1^- and L_2 matrices are defined as$\\begin{aligned} L1^{-} &\\equiv \\frac{1}{\\Delta}\\begin{pmatrix} 1 - (1 + \\xi \\Delta) &0&0&\\dots&0&0&0\\\n-1&1&0&\\dots&0&0&0\\\n\\vdots&\\vdots&\\vdots&\\ddots&\\vdots&\\vdots&\\vdots\\\n0&0&0&\\dots&-1&1&0\\\n0&0&0&\\cdots&0&-1&1 \\end{pmatrix}{P\\times P}\\\nL2 &\\equiv \\frac{1}{\\Delta^2}\\begin{pmatrix} -2 + (1 + \\xi\\Delta) &1&0&\\dots&0&0&0\\\n1&-2&1&\\dots&0&0&0\\\n\\vdots&\\vdots&\\vdots&\\ddots&\\vdots&\\vdots&\\vdots\\\n0&0&0&\\dots&1&-2&1\\\n0&0&0&\\cdots&0&1&-2 + (1- \\xi\\Delta) \\end{pmatrix}{P\\times P} \\end{aligned} $"
},

{
    "location": "#Irregular-grids-1",
    "page": "﻿Installation",
    "title": "Irregular grids",
    "category": "section",
    "text": "Using the backward first-order difference, the first boundary condition (1) can be alternatively represented as $\\begin{aligned} \\dfrac{v(0) - v(-\\Delta{1, -})}{\\Delta{1, -}} &= - \\xi v(0) \\end{aligned} $Note that we have assumed that Delta_1- = Delta_1+ and Delta_P+ = Delta_P- for the ghost notes. Using discretized central differences of second orders, the first boundary condition (1) yields$\\begin{aligned} &\\dfrac{\\Delta{1,-} v( \\Delta{1,+}) - (\\Delta{1,+} + \\Delta{1,-}) v( 0 ) + \\Delta{1,+}  v( - \\Delta{1,-})}{\\frac{1}{2}(\\Delta{1,+} + \\Delta{1,-}) \\Delta{1,+} \\Delta{1,-} } \\\n&= \\dfrac{v (\\Delta{1, +}) - 2 v(0) + v(-\\Delta{1, +})}{\\Delta{1, +}^2} \\ &= \\dfrac{v(\\Delta{1, +}) - v(0)}{\\Delta{1, +}^2} - \\dfrac{1}{\\Delta{1, +}}\\dfrac{v (0) - v(-\\Delta{1, +}) }{\\Delta{1, +}}  \\\n&= \\dfrac{v(\\Delta{1, +}) - v(0)}{\\Delta{1, +}^2} + \\dfrac{1}{\\Delta{i,+}} \\xi v(0)  \\  &= \\dfrac{1}{\\Delta{1, +}^2}  (- 1 + \\Delta{1, +} \\xi) v(0)  + \\dfrac{1}{\\Delta{1, +}^2}  v(\\Delta_{1, +})   \\end{aligned} $Similarly, for the secondary boundary condition (2), we have $\\begin{aligned} &\\dfrac{\\Delta{P,-} v( \\bar{z} + \\Delta{P,+}) - (\\Delta{P,+} + \\Delta{P,-}) v(\\bar{z} ) + \\Delta{P,+}  v( \\bar{z} - \\Delta{P,-})}{\\frac{1}{2}(\\Delta{P,+} + \\Delta{P,-}) \\Delta{P,+} \\Delta{P,-} } \\\n&=\\dfrac{v (\\bar{z} + \\Delta{P,-}) - 2 v(\\bar{z} ) + v(\\bar{z} -\\Delta{P,-})}{\\Delta{P,-}^2} \\\n&=   \\dfrac{v(\\bar{z} - \\Delta{P,-}) - v(\\bar{z})}{\\Delta{P,-}^2} + \\dfrac{1}{\\Delta{P,-}}\\dfrac{ v(\\bar{z}+\\Delta{P,-}) - v (\\bar{z}) }{\\Delta{P,-}}  \\\n&= \\dfrac{v(\\bar{z} - \\Delta{P,-}) - v(\\bar{z})}{\\Delta{P,-}^2}  - \\dfrac{1}{\\Delta{P,-}} \\xi v(\\bar{z})  \\  &= \\dfrac{1}{\\Delta{P,-}^2}  (- 1 - \\Delta{P,-} \\xi) v(\\bar{z})  + \\dfrac{1}{\\Delta{P,-}^2}  v(\\bar{z} - \\Delta_{P,-})   \\end{aligned} $Thushe corresponding L_1^- and L_2 matrices are defined as$\\begin{aligned} L1^{-} &\\equiv \\begin{pmatrix} \\Delta^{-1}{1,-} [1 - (1 + \\xi \\Delta^{-1}{1,-})] &0&0&\\dots&0&0&0\\\n-\\Delta{2,-}^{-1}&\\Delta{2,-}^{-1}&0&\\dots&0&0&0\\\n\\vdots&\\vdots&\\vdots&\\ddots&\\vdots&\\vdots&\\vdots\\\n0&0&0&\\dots&-\\Delta{P-1,-}^{-1}&\\Delta{P-1,-}^{-1}&0\\\n0&0&0&\\cdots&0&-\\Delta{P,-}^{-1}&\\Delta{P,-}^{-1} \\end{pmatrix}{P\\times P} \\\nL2 &\\equiv \\begin{pmatrix} \\Delta{1,+}^{-2}[-2 + (1+\\xi \\Delta{1,+})] &\\Delta{1,+}^{-2}&0&\\cdots&0&0&0 \\\n\\vdots&\\ddots&\\ddots&\\ddots&\\ddots&\\vdots&\\vdots\\\n0&\\cdots&2(\\Delta{i,+}+\\Delta{i,-})^{-1} \\Delta{i,-}^{-1} &-2\\Delta{i,-}^{-1} \\Delta{i,+}^{-1}  & 2 (\\Delta{i,+}+\\Delta{i,-})^{-1} \\Delta{i,+}^{-1}&\\cdots&0 \\\n\\vdots&\\vdots&\\vdots&\\ddots&\\ddots&\\ddots&\\vdots\\\n0&0&0&\\cdots&0&\\Delta{P,-}^{-2}&\\Delta{P,-}^{-2} [-2 + (1- \\xi\\Delta{P,-})] \\end{pmatrix}{P\\times P} \\end{aligned} $"
},

]}
