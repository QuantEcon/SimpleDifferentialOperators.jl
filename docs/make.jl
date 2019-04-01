using Documenter, SimpleDifferentialOperators

# Compile the online documentation.
makedocs(sitename = "SimpleDifferentialOperators.jl",
	pages = [
        "index.md",
        "notebooks.md",
        "api.md"
    ],
	doctest = :fix)

# Deploy to git
deploydocs(
    repo = "github.com/QuantEcon/SimpleDifferentialOperators.jl.git",
)
