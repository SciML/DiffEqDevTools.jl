using Documenter, DiffEqDevTools

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "DiffEqDevTools.jl",
         authors = "Chris Rackauckas",
         modules = [DiffEqDevTools],
         clean = true,
         doctest = false,
         format = Documenter.HTML(assets = ["assets/favicon.ico"],
                                  canonical = "https://docs.sciml.ai/DiffEqDevTools/stable/"),
         pages = pages)

deploydocs(repo = "github.com/SciML/DiffEqDevTools.jl.git"; push_preview = true)
