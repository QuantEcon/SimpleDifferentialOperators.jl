using Documenter, SimpleDifferentialOperators

# Compile the online documentation.
makedocs(sitename = "SimpleDifferentialOperators.jl",
	pages = [
        "index.md",
        "formula.md",
        "api.md",
    ])

# Deploy to git
deploydocs(
    repo = "github.com/QuantEcon/SimpleDifferentialOperators.jl.git",
)
