using Documenter, SimpleDifferentialOperators

# Compile the raw documentation.
makedocs(sitename = "SimpleDifferentialOperators.jl")

# Deploy to git
deploydocs(
    repo = "github.com/QuantEcon/SimpleDifferentialOperators.jl.git",
)
