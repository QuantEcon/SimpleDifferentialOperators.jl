using Weave
fileset = readdir(@__DIR__)
for file in fileset
    occursin(".jmd", file) && Weave.notebook(file, @__DIR__) && Weave.weave(file, out_path = @__DIR__, doctype = "md2html")
    println("Weaved $file")
end
