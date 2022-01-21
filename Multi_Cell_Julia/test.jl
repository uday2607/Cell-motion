area_clib = joinpath(@__DIR__, "./C_area_codes/libmyclib.so")

function c_ell_area(a1, b1, h1,
    k1, t1, a2, b2,
    h2, k2, t2)

    area = ccall((:ee_ovlap_area, area_clib), Float64, (Float64, Float64, 
            Float64, Float64, Float64, Float64, Float64, Float64, 
            Float64, Float64), a1, b1, h1, k1, t1, a2, b2, h2, k2, t2)

    return area
end    