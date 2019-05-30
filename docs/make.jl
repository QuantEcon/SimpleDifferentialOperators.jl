using Documenter, SimpleDifferentialOperators

# figure out what to paste instead of TAG_GOES_HERE 
if "TRAVIS_TAG" ∈ keys(ENV)
    replacement = get(ENV, "TRAVIS_TAG")
else 
    replacement = "dev"
end 

fileset = readdir(@__DIR__)
for file in fileset
    if occursin(".jmd", file) 
        txt = read(file, String)
        open(file, "w") do f
           write(f, replace(txt, "TAG_GOES_HERE" => replacement))
        end     
    end
end

# Compile the online documentation.
makedocs(sitename = "SimpleDifferentialOperators.jl",
	pages = [
        "Home" => "index.md", 
        "Gallery" => [
            "Examples" => "examples.md",
            "Notebooks" => "notebooks.md"
        ],
        "api.md"
    ],
	doctest = :fix)

# Deploy to git
deploydocs(
    repo = "github.com/QuantEcon/SimpleDifferentialOperators.jl.git",
)

