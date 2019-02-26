using Documenter, SimpleDifferentialOperators, Weave

# Generate stuff from Weave
srcfile = "example.jmd" # example script from weave docs
weave(srcfile, doctype = "md2pdf", out_path = :pwd) # turn it into PDF
convert_doc(srcfile, "example.ipynb") # try it as notebook?

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
