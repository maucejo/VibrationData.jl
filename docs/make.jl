using VibrationData
using Documenter

makedocs(;
    # modules=[VibrationData],
    authors="Mathieu Aucejo",
    sitename="VibrationData.jl",
    format=Documenter.HTML(;
        canonical="https://maucejo.github.io/VibrationData.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/maucejo/VibrationData.jl",
    devbranch="master",
)
