include("corner-algorithms.jl")

"""
```
corners = imcorner(img; [method])
corners = imcorner(img, threshold; [method])
```
Performs corner detection using one of the following methods -
    1. harris
    2. shi_tomasi
    3. kitchen_rosenfeld
The parameters of the individual methods are described in their documentation. The
maxima values of the resultant responses are taken as corners. If a threshold is
specified, the values of the responses are thresholded to give the corner pixels.
If `threshold` is a `Percentile` then its type will be preserved.
"""
function imcorner(img::AbstractArray; method::Function = harris, args...)
    responses = method(img; args...)
    corners = similar(img, Bool)
    fill!(corners, false)
    maxima = map(CartesianIndex, findlocalmaxima(responses))
    for m in maxima corners[m] = true end
    corners
end

function imcorner(img::AbstractArray, threshold; method::Function = harris, args...)
    responses = method(img; args...)
    map(i -> i > threshold, responses)
end

function imcorner(img::AbstractArray, thresholdp::Percentile; method::Function = harris, args...)
    responses = method(img; args...)
    threshold = StatsBase.percentile(vec(responses), thresholdp.p)
    map(i -> i > threshold, responses)
end

"""
```
corners = imcorner_subpixel(img; [method])
         -> Vector{HomogeneousPoint{Float64,3}}
corners = imcorner_subpixel(img, threshold, percentile; [method])
         -> Vector{HomogeneousPoint{Float64,3}}
```
Same as [`imcorner`](@ref), but estimates corners to sub-pixel precision.
Sub-pixel precision is achieved by interpolating the corner response values using
the 4-connected neighbourhood of a maximum response value.
See [`corner2subpixel`](@ref) for more details of the interpolation scheme.
"""
function imcorner_subpixel(img::AbstractArray; method::Function = harris, args...)
    responses = method(img; args...)
    corner_indicator = similar(img, Bool)
    fill!(corner_indicator, false)
    maxima = map(CartesianIndex, findlocalmaxima(responses))
    for m in maxima corner_indicator[m] = true end
    corners = corner2subpixel(responses,corner_indicator)
end

function imcorner_subpixel(img::AbstractArray, threshold; method::Function = harris, args...)
    responses = method(img; args...)
    corner_indicator = map(i -> i > threshold, responses)
    corners = corner2subpixel(responses,corner_indicator)
end

function imcorner_subpixel(img::AbstractArray, thresholdp::Percentile; method::Function = harris, args...)
    responses = method(img; args...)
    threshold = StatsBase.percentile(vec(responses), thresholdp.p)
    corner_indicator = map(i -> i > threshold, responses)
    corners = corner2subpixel(responses,corner_indicator)
end


"""
```
corners = corner2subpixel(responses::AbstractMatrix,corner_indicator::AbstractMatrix{Bool})
        -> Vector{HomogeneousPoint{Float64,3}}
```
Refines integer corner coordinates to sub-pixel precision.
The function takes as input a matrix representing corner responses
and a boolean indicator matrix denoting the integer coordinates of a corner
in the image. The output is a vector of type [`HomogeneousPoint`](@ref)
storing the sub-pixel coordinates of the corners.
The algorithm computes a correction factor which is added to the original
integer coordinates. In particular, a univariate quadratic polynomial is fit
separately to the ``x``-coordinates and ``y``-coordinates of a corner and its immediate
east/west, and north/south neighbours. The fit is achieved using a local
coordinate system for each corner, where the origin of the coordinate system is
a given corner, and its immediate neighbours are assigned coordinates of  minus
one and plus one.
The corner and its two neighbours form a system of three equations. For example,
let  ``x_1 = -1``,  ``x_2 = 0`` and  ``x_3 = 1`` denote the local ``x`` coordinates
of the west, center and east pixels and let the vector ``\\mathbf{b} = [r_1, r_2, r_3]``
denote the corresponding corner response values. With
```math
    \\mathbf{A} =
        \\begin{bmatrix}
            x_1^2 & x_1  & 1  \\\\
            x_2^2 & x_2  & 1 \\\\
            x_3^2 & x_3  & 1 \\\\
        \\end{bmatrix},
```
the coefficients of the quadratic polynomial can be found by solving the
system of equations ``\\mathbf{b} = \\mathbf{A}\\mathbf{x}``.
The result is given by ``x = \\mathbf{A}^{-1}\\mathbf{b}``.
The vertex of the quadratic polynomial yields a sub-pixel estimate of the
true corner position. For example, for a univariate quadratic polynomial
``px^2 + qx + r``, the ``x``-coordinate of the vertex is ``\\frac{-q}{2p}``.
Hence, the refined sub-pixel coordinate is equal to:
 ``c +  \\frac{-q}{2p}``, where ``c`` is the integer coordinate.
!!! note
    Corners on the boundary of the image are not refined to sub-pixel precision.
"""
function corner2subpixel(responses::AbstractMatrix, corner_indicator::AbstractMatrix{Bool})
    row_range, col_range = axes(corner_indicator)
    idxs = findall(!iszero, corner_indicator) # findnz
    row, col = (getindex.(idxs,1),getindex.(idxs,2))
    ncorners = length(row)
    corners = fill(HomogeneousPoint((0.0,0.0,0.0)),ncorners)
    invA = @SMatrix [0.5 -1.0 0.5; -0.5 0.0 0.5; 0.0 1.0 -0.0]
    for k = 1:ncorners
        # Corners on the perimeter of the image will not be interpolated.
        if  (row[k] == first(row_range) || row[k] == last(row_range) ||
             col[k] == first(col_range) || col[k] == last(col_range))
            y = convert(Float64,row[k])
            x = convert(Float64,col[k])
            corners[k] = HomogeneousPoint((x,y,1.0))
        else
            center, north, south, east, west =
                                unsafe_neighbourhood_4(responses,row[k],col[k])
            # Solve for the coefficients of the quadratic equation.
            a, b, c = invA* @SVector [west, center, east]
            p, q, r = invA* @SVector [north, center, south]
            # Solve for the first coordinate of the vertex.
            u = -b/(2.0a)
            v = -q/(2.0p)
            corners[k] = HomogeneousPoint((col[k]+u,row[k]+v,1.0))
        end
    end
    return corners
end
