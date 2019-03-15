using Weave
fileset = readdir(@__DIR__)
for file in fileset
    if occursin(".jmd", file)
        Weave.notebook(file)
        Weave.weave(file, out_path = :pwd, doctype = "md2html", template = "our_template.tpl")
        println("Weaved $file")
    end
end
