push!(LOAD_PATH,"../src/")
using ImageCorners
using Documenter

DocMeta.setdocmeta!(ImageCorners, :DocTestSetup, :(using ImageCorners); recursive=true)

makedocs(;
    modules=[ImageCorners],
    sitename="ImageCorners.jl",
)

deploydocs(;
    repo="github.com/JuliaImages/ImageCorners.jl",
)