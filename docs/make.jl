using SpectralBoundaryIntegralMethod
using Documenter

DocMeta.setdocmeta!(SpectralBoundaryIntegralMethod, :DocTestSetup, :(using SpectralBoundaryIntegralMethod); recursive=true)

makedocs(;
    modules=[SpectralBoundaryIntegralMethod],
    authors="Ali Gurbuz",
    repo="https://github.com/aligurbu/SpectralBoundaryIntegralMethod.jl/blob/{commit}{path}#{line}",
    sitename="SpectralBoundaryIntegralMethod.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
