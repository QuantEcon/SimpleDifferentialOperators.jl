using Weave
fileset = readdir(@__DIR__)
for file in fileset
    occursin(".jmd", file) && Weave.notebook(file, @__DIR__)
    occursin(".jmd", file) && Weave.weave(file, out_path = :pwd, doctype = "md2html", template = "our_template.tpl")
    println("Weaved $file")
end
