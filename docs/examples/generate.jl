using Weave
fileset = readdir(@__DIR__)
for file in fileset
    occursin(".jmd", file) && Weave.notebook(file, @__DIR__)
    println("Weaved $file")
end
