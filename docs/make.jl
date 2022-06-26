using NonstationaryProcesses
using Documenter

DocMeta.setdocmeta!(NonstationaryProcesses, :DocTestSetup, :(using NonstationaryProcesses); recursive=true)

makedocs(;
    modules=[NonstationaryProcesses],
    authors="brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
    repo="https://github.com/brendanjohnharris/NonstationaryProcesses.jl/blob/{commit}{path}#{line}",
    sitename="NonstationaryProcesses.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
