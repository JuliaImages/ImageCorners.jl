module ImageCorners

using ImageCore
using ImageCore: NumberLike
using ImageFiltering
using PrecompileTools
using StaticArrays
using StatsBase

include("utils.jl")
include("cornerapi.jl")
include("precompile.jl")

export
    imcorner,
    imcorner_subpixel,
    corner2subpixel,
    harris,
    shi_tomasi,
    kitchen_rosenfeld,
    moravec,
    fastcorners,
    meancovs,
    gammacovs,
   

    Percentile,
    HomogeneousPoint

end # module ImageCorners
