var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "To install, run] add SimpleDifferentialOperatorsNote that this requires Julia 1.0 or later.For the complete derivations of the objects we return, see the PDF."
},

{
    "location": "#Formulas-1",
    "page": "Installation",
    "title": "Formulas",
    "category": "section",
    "text": "Under P-length of grids on v with end points underlinez  overlinez with reflecting barrier conditions ofbeginalign\nunderlinexi v(underlinez) + nabla v(underlinez ) = 0\noverlinexi v(overlinez) + nabla v(overlinez) = 0\nendalignWe use the following discretization schemes:L_1^- equiv frac1Deltabeginpmatrix\n1 - (1 + underlinexi Delta) 00dots000\n-110dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Ptimes PL_1^+ equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots0-11\n000cdots00-1+(1-overlinexi Delta)\nendpmatrix_Ptimes PlabeleqL-1-plus-regular L_2 equiv frac1Delta^2beginpmatrix\n-2 + (1 + underlinexi Delta) 10dots000\n1-21dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots1-21\n000cdots01-2 + (1- overlinexi Delta)\nendpmatrix_Ptimes PlabeleqL-2-regularwhich represent the backward first order, foward first order, and central second order differential operators respectively.Derivation, including formula for irregular grids, can be found here."
},

{
    "location": "#SimpleDifferentialOperators.diffusionoperators-Tuple{Any,BoundaryCondition,BoundaryCondition}",
    "page": "Installation",
    "title": "SimpleDifferentialOperators.diffusionoperators",
    "category": "method",
    "text": "diffusionoperators(x, BC1::BoundaryCondition, BC2::BoundaryCondition)\n\nReturn the diffusion operators (L_1_minus, L_1_plus, L_2) w.r.t the supplied grid and BCs.\n\nx is a grid (either an AbstractRange, in which we use specialized uniform grid code, or an AbstractArray). The first BC binds at the lower end of the grid (i.e., x[1]), and the latter at the high end. The BCs are either a Reflecting(), or \"Dirichlet\" boundary condition v\'(x) = 0, or Mixed(x::T) where {T <: Real}, corresponding to \"Robin\" boundary conditions.\n\n\n\n\n\n"
},

{
    "location": "#Usage-1",
    "page": "Installation",
    "title": "Usage",
    "category": "section",
    "text": "This package provides discretized differential operators of first order and second order under reflecting boundary conditions.Modules = [SimpleDifferentialOperators]\nOrder   = [:function, :type]"
},

]}
