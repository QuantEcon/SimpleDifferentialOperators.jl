using Documenter, SimpleDifferentialOperators, Weave

# Generate stuff from Weave
srcfile = joinpath(pwd(), "docs", "examples", "example.jmd") # example script from weave docs
Weave.weave(srcfile, doctype = "md2html", out_path = joinpath(pwd(), "docs", "build", "examples")) # turn it into PDF
convert_doc(srcfile, "docs/build/examples/example.ipynb") # try it as notebook?

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
