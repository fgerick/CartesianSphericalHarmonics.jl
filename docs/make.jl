using CartesianSphericalHarmonics
using Documenter

makedocs(;
    modules=[CartesianSphericalHarmonics],
    authors="Felix Gerick",
    repo="https://github.com/fgerick/CartesianSphericalHarmonics.jl/blob/{commit}{path}#L{line}",
    sitename="CartesianSphericalHarmonics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fgerick.github.io/CartesianSphericalHarmonics.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fgerick/CartesianSphericalHarmonics.jl",
)
