using PrecompileTools

@info "Precompiling ImageCorners"

@setup_workload begin

    imgs_list = Any[
        rand(Gray{ColorTypes.N0f8}, 32, 32),
        rand(RGB{ColorTypes.N0f8}, 32, 32),
        rand(Gray{ColorTypes.Float64}, 32, 32),
        rand(RGB{Float64}, 32, 32),
        rand(ColorTypes.N0f8, 32, 32),
        rand(ColorTypes.Float64, 32, 32)
    ]

    @compile_workload begin
        for img in imgs_list
            imcorner(img, method = harris)
            imcorner(img, method = shi_tomasi)
            imcorner(img, method = kitchen_rosenfeld)

            imcorner(img, 0.001, method = harris)
            imcorner(img, 0.001, method = shi_tomasi)
            imcorner(img, 0.001, method = kitchen_rosenfeld)

            imcorner(img, Percentile(95), method = harris)
            imcorner(img, Percentile(95), method = shi_tomasi)
            imcorner(img, Percentile(95), method = kitchen_rosenfeld)

            imcorner_subpixel(img, method = harris)
            imcorner_subpixel(img, method = shi_tomasi)
            imcorner_subpixel(img, method = kitchen_rosenfeld)
           
            imcorner_subpixel(img, 0.001, method = harris)
            imcorner_subpixel(img, 0.001, method = shi_tomasi)
            imcorner_subpixel(img, 0.001, method = kitchen_rosenfeld)

            imcorner_subpixel(img, Percentile(95), method = harris)
            imcorner_subpixel(img, Percentile(95), method = shi_tomasi)
            imcorner_subpixel(img, Percentile(95), method = kitchen_rosenfeld)

            harris(img)
            shi_tomasi(img)
            kitchen_rosenfeld(img)

            fastcorners(img)
            fastcorners(img, 10)
            fastcorners(img, 10, 1)
        end

        Ac = zeros(41,41)
        I = fill(false,41,41)
        corner2subpixel(Ac, I)
    end
end