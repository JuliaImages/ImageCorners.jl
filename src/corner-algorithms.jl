"""
```
harris_response = harris(img; [k], [border], [weights])
```
Performs Harris corner detection. The covariances can be taken using either a mean
weighted filter or a gamma kernel.
"""
function harris(img::AbstractArray; k::Float64 = 0.04, args...)
    cov_xx, cov_xy, cov_yy = gradcovs(img, args...)
    corner = map((xx, yy, xy) -> xx * yy - xy ^ 2 - k * (xx + yy) ^ 2, cov_xx, cov_yy, cov_xy)
    corner
end

"""
```
shi_tomasi_response = shi_tomasi(img; [border], [weights])
```
Performs Shi Tomasi corner detection. The covariances can be taken using either a mean
weighted filter or a gamma kernel.
"""
function shi_tomasi(img::AbstractArray; border::AbstractString = "replicate", args...)
    cov_xx, cov_xy, cov_yy = gradcovs(img, border; args...)
    corner = map((xx, yy, xy) -> ((xx + yy) - (sqrt((xx - yy) ^ 2 + 4 * xy ^ 2))) / 2, cov_xx, cov_yy, cov_xy)
    corner
end

"""
```
kitchen_rosenfeld_response = kitchen_rosenfeld(img; [border])
```
Performs Kitchen Rosenfeld corner detection. The covariances can be taken using either a mean
weighted filter or a gamma kernel.
"""
function kitchen_rosenfeld(img::AbstractArray; border::AbstractString = "replicate")
    meth = KernelFactors.sobel
    (grad_x, grad_y) = imgradients(img, meth, border)
    (grad_xx, grad_xy) = imgradients(grad_x, meth, border)
    (grad_yx, grad_yy) = imgradients(grad_y, meth, border)
    map(kr, grad_x, grad_y, grad_xx, grad_xy, grad_yy)
end

function kr(x::T, y::T, xx::T, xy::T, yy::T) where T<:Real
    num = xx*y*y + yy*x*x - 2*xy*x*y
    denom = x*x + y*y
    ifelse(denom == 0, zero(num)/one(denom), -num/denom)
end

function kr(x::Real, y::Real, xx::Real, xy::Real, yy::Real)
    xp, yp, xxp, xyp, yyp = promote(x, y, xx, xy, yy)
    kr(xp, yp, xxp, xyp, yyp)
end

kr(x::NumberLike, y::NumberLike, xx::NumberLike, xy::NumberLike, yy::NumberLike) =
    kr(gray(x), gray(y), gray(xx), gray(xy), gray(yy))

function kr(x::AbstractRGB, y::AbstractRGB, xx::AbstractRGB, xy::AbstractRGB, yy::AbstractRGB)
    krrgb = RGB(kr(red(x), red(y), red(xx), red(xy), red(yy)),
                kr(green(x), green(y), green(xx), green(xy), green(yy)),
                kr(blue(x),  blue(y),  blue(xx),  blue(xy),  blue(yy)))
    gray(convert(Gray, krrgb))
end

"""
    fastcorners(img, n, threshold) -> corners
Performs FAST Corner Detection. `n` is the number of contiguous pixels
which need to be greater (lesser) than intensity + threshold (intensity - threshold)
for a pixel to be marked as a corner. The default value for n is 12.
"""
function fastcorners(img::AbstractArray{T}, n::Int = 12, threshold::Float64 = 0.15) where T
    img_padded = padarray(img, Fill(0, (3,3)))
    corner = falses(size(img))
    R = CartesianIndices(size(img))
    idx = map(CartesianIndex{2}, [(0, 3), (1, 3), (2, 2), (3, 1), (3, 0), (3, -1), (2, -2), (1, -3),
            (0, -3), (-1, -3), (-2, -2), (-3, -1), (-3, 0), (-3, 1), (-2, 2), (-1, 3)])

    idxidx = [1, 5, 9, 13]
    for I in R
        bright_threshold = img_padded[I] + threshold
        dark_threshold = img_padded[I] - threshold
        if n >= 12
            sum_bright = 0
            sum_dark = 0
            for k in idxidx
                pixel = img_padded[I + idx[k]]
                if pixel > bright_threshold
                    sum_bright += 1
                elseif pixel < dark_threshold
                    sum_dark += 1
                end
            end
            if sum_bright < 3 && sum_dark < 3
                continue
            end
        end
        consecutive_bright = 0
        consecutive_dark = 0

        for i in 1:15 + n
            k = mod1(i, 16)
            pixel = img_padded[I + idx[k]]
            if pixel > bright_threshold
                consecutive_dark = 0
                consecutive_bright += 1
            elseif pixel < dark_threshold
                consecutive_bright = 0
                consecutive_dark += 1
            end

            if consecutive_dark == n || consecutive_bright == n
                corner[I] = true
                break
            end
        end
    end
    corner
end

function gradcovs(img::AbstractArray, border::AbstractString = "replicate"; weights::Function = meancovs, args...)
    (grad_x, grad_y) = imgradients(img, KernelFactors.sobel, border)

    cov_xx = dotc.(grad_x, grad_x)
    cov_xy = dotc.(grad_x, grad_y)
    cov_yy = dotc.(grad_y, grad_y)

    weights(cov_xx, cov_xy, cov_yy, args...)
end

function meancovs(cov_xx, cov_xy, cov_yy, blockSize::Int = 3)

    box_filter_kernel = centered((1 / (blockSize * blockSize)) * ones(blockSize, blockSize))

    filt_cov_xx = imfilter(cov_xx, box_filter_kernel)
    filt_cov_xy = imfilter(cov_xy, box_filter_kernel)
    filt_cov_yy = imfilter(cov_yy, box_filter_kernel)

    filt_cov_xx, filt_cov_xy, filt_cov_yy
end

function gammacovs(cov_xx, cov_xy, cov_yy, gamma::Float64 = 1.4)
    kernel = KernelFactors.gaussian((gamma, gamma))

    filt_cov_xx = imfilter(cov_xx, kernel)
    filt_cov_xy = imfilter(cov_xy, kernel)
    filt_cov_yy = imfilter(cov_yy, kernel)

    filt_cov_xx, filt_cov_xy, filt_cov_yy
end


"""
```
moravec_response = moravec(img; [threshold], [window_size])
```
Compute the corner response values for each pixel in the input image using the Moravec corner detection algorithm.

    # Arguments
    - `img::AbstractArray`: The input image.
    - `window_size::Int`: Size of the window for calculating the corner response.
    
    # Keyword Arguments
    - `args...`: Additional arguments passed to the `gradcovs` function.
    # Resources
    - https://vincmazet.github.io/bip/detection/corners.html
    - http://www0.cs.ucl.ac.uk/staff/g.brostow/classes/IP2008/L7_CornerDetection.pdf

"""
function moravec(img::AbstractArray, args...; threshold::Float64 = 10000.0, window_size::Int = 3)
    gradient_x, gradient_y = gradcovs(img, args...)

    corners = []
    for y in 1:size(img, 1) - window_size + 1
        for x in 1:size(img, 2) - window_size + 1
            min_sum_diff = Inf
            for dy in -1:1
                for dx in -1:1
                    sum_diff = sum((gradient_x[y:y+window_size-1, x:x+window_size-1] .- gradient_x[y+dy:y+dy+window_size-1, x+dx:x+dx+window_size-1]).^2)
                    min_sum_diff = min(min_sum_diff, sum_diff)
                end
            end
            img[y,x] = min_sum_diff
        end
    end

end
