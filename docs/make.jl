using Documenter, alfa

makedocs(;
    modules=[alfa],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/NilsKintscher/alfa.jl/blob/{commit}{path}#L{line}",
    sitename="alfa.jl",
    authors="Nils Kintscher, Karsten Kahl",
    assets=String[],
)

deploydocs(;
    repo="github.com/NilsKintscher/alfa.jl",
)
