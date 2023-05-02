# ImageCorners.jl

ImageCorners.jl provides image corner related algorithms in Julia. 

### Supported Algorithms
- Harris 
- Shi-Tomasi 
- Kitchen and Rosenfeld 
- FAST Corners

These Algorithms can be accessed using the following methods:
```julia
julia> using ImageCorners, TestImages

julia> img = testimage("mandrill");

# corner detection using harris method
julia> corners = imcorner(img; method = harris);

# threshold can be specified for the thresolding the corner pixels
julia> corners = imcorner(img, 0.001; method = harris);

julia> corners = imcorner(img, Percentile(95); method = harris);

# for corner detection to subpixel precision imgcorner_subpixel can be used
julia> corners = imcorner_subpixel(img; method = harris);
```
