using Documenter, SimpleDifferentialOperators, Weave, IJulia

# Generated files
EXAMPLE = joinpath(@__DIR__, "..", "examples", "example.jmd")
OUTPUT = joinpath(@__DIR__, "src", "generated", "example.ipynb")
convert_doc(EXAMPLE, OUTPUT)
Weave.notebook(OUTPUT)

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
