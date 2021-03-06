using OptimalBids
using Documenter

DocMeta.setdocmeta!(OptimalBids, :DocTestSetup, :(using OptimalBids); recursive=true)

makedocs(;
    modules=[OptimalBids],
    authors="Andrew Rosemberg",
    repo="https://github.com/andrewrosemberg/OptimalBids.jl/blob/{commit}{path}#{line}",
    sitename="OptimalBids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://andrewrosemberg.github.io/OptimalBids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Nonconvex.jl" => "examples/Nonconvex.md",
            "Optim.jl" => "examples/Optim.md",
            "AbstractGPs.jl" => "examples/GaussianProcesses.md",
        ],
        "API" => "api.md",
    ],
    checkdocs=:exports,
    strict=false,
)

deploydocs(; repo="github.com/andrewrosemberg/OptimalBids.jl", devbranch="main")
