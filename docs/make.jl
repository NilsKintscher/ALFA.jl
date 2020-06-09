using Documenter, ALFA, ALFA.gallery

makedocs(;
    modules=[ALFA, ALFA.gallery],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Examples" =>  map(
                s -> "examples/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/examples")))
            ),
        "Library" =>  map(
                s -> "internals/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/internals")))
            ),
    ],
    repo="https://github.com/NilsKintscher/ALFA.jl/blob/{commit}{path}#L{line}",
    sitename="ALFA.jl",
    authors="Nils Kintscher, Karsten Kahl",
)

deploydocs(;
    repo="github.com/NilsKintscher/ALFA.jl",
)
