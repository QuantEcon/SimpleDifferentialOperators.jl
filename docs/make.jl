using Pkg
pkg"instantiate"
using Documenter, SimpleDifferentialOperators, Weave

# Generated files
EXAMPLE = joinpath(@__DIR__, "..", "examples", "example.jmd")
OUTPUT = joinpath(@__DIR__, "src", "generated", "example.ipynb")
Weave.notebook(EXAMPLE, OUTPUT)

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
