using Documenter, SimpleDifferentialOperators

# Compile the raw documentation.
makedocs(sitename = "SimpleDifferentialOperators.jl")

# MkDocs
deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/QuantEcon/SimpleDifferentialOperators.jl.git",
    target = "site",
    make = () -> run(`mkdocs build`)
)
