using Documenter
push!(LOAD_PATH, "../src/")
using JLACovarianceMatrix

DocMeta.setdocmeta!(JLACovarianceMatrix, :DocTestSetup, :(using JLACovarianceMatrix); recursive=true)

makedocs(
    sitename="JLACovarianceMatrix Documentation",
    modules = [JLACovarianceMatrix],
    pages = [
        "JLACovarianceMatrix" => "index.md",
    ],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
    )
)

deploydocs(
    repo = "github.com/OmegaLambda1998/JLACovarianceMatrix.jl.git"
)
