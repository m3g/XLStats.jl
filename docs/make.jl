import Pkg
Pkg.add("Documenter")
using Documenter
using XLStats
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[XLStats],
    sitename="XLStats.jl",
    pages = [
        "Home" => "index.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/XLStats.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
