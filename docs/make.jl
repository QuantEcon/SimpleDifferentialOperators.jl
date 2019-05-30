using Documenter, SimpleDifferentialOperators

# figure out what to paste instead of TAG_GOES_HERE 
replacement = get(ENV, "TRAVIS_TAG", "dev")

tmp = @__DIR__
fileset = readdir(tmp * "/src/")
for file in fileset
    if occursin(".md", file) 
        txt = read(tmp * "/src/" * file, String)
        open(tmp * "/src/" * file, "w") do f
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

