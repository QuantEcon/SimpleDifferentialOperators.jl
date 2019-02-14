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
    "location": "#Usage-1",
    "page": "Installation",
    "title": "Usage",
    "category": "section",
    "text": "This package provides discretized differential operators of first order and second order under reflecting boundary conditions."
},

{
    "location": "#Formula-1",
    "page": "Installation",
    "title": "Formula",
    "category": "section",
    "text": "We use the following discretization schemes on P-length of grids on v with end points underlinez  barz under reflecting barrier conditions ofbeginalign\nxi v(underlinez) + Dzv(underlinez ) = 0\nxi v(barz) + Dzv(barz) = 0\nendalignL_1^- equiv frac1Deltabeginpmatrix\n1 - (1 + xi Delta) 00dots000\n-110dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots-110\n000cdots0-11\nendpmatrix_Ptimes PL_1^+ equiv frac1Deltabeginpmatrix\n-110dots000\n0-11dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots0-11\n000cdots00-1+(1-xi Delta)\nendpmatrix_Ptimes PlabeleqL-1-plus-regular L_2 equiv frac1Delta^2beginpmatrix\n-2 + (1 + xiDelta) 10dots000\n1-21dots000\nvdotsvdotsvdotsddotsvdotsvdotsvdots\n000dots1-21\n000cdots01-2 + (1- xiDelta)\nendpmatrix_Ptimes PlabeleqL-2-regularwhich represent the backward first order, foward first order, and central second order differential operators respectively.Derivation, including formula for irregular grids, can be found here."
},

]}
