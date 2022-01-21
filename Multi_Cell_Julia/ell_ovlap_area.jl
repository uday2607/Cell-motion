using PolygonArea
using StaticArrays

const Point = SVector{2}

# Ellipses 
function ellipse(x0, y0, a, b, phi, nb_vertices)
    vertices = [Point(x0 + a*cos(θ)*cos(phi) - b*sin(θ)*sin(phi), 
                         y0 + a*cos(θ)*sin(phi) + b*sin(θ)*cos(phi)) 
                         for θ in reverse(
                         LinRange(0.0, 2π, nb_vertices+1)[1:nb_vertices])]
    # Reversed, such that the inner of the domain is on the right of the interface
    return ConvexPolygon(vertices)
end

function ellipse_ellipse_overlap(a1, b1, h1,
    k1, t1, a2, b2,
    h2, k2, t2)

    ell1 = ellipse(h1, k1, a1, b1, t1, 45)
    ell2 = ellipse(h2, k2, a2, b2, t2, 45)

    Area = area(ell1 ∩ ell2)

    return Area
end    

function ellipse_ellipse_overlap_e(ell1, ell2)

    Area = area(ell1 ∩ ell2)

    return Area
end 