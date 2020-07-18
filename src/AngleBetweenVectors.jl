module AngleBetweenVectors

export arc

import Base: angle

import LinearAlgebra: norm


@inline unitize(p) = p ./ norm(p)

function arc(point1::NTuple{2,T}, point2::NTuple{2,T}) where {T<:Real}
    return atan(point1[2], point1[1]) - atan(point2[2], point2[1])
end

function arc(pt1::NTuple{3,T}, pt2::NTuple{3,T}) where {T<:Real}
    n1 = ab_minus_cd(pt1[2],pt2[3],pt1[3],pt2[2])
    n2 = ab_minus_cd(pt1[3],pt2[1],pt1[1],pt2[3])
    n3 = ab_minus_cd(pt1[1],pt2[2],pt1[2],pt2[1])
    h = hypot(n1,n2,n3)
    n1 = n1 / h; n2 = n2 / h; n3 = n3 / h
    h = hypot(n1,n2,n3)
    return atan(h)
end


# from Kahan
function ad_minus_bc(a,b,c,d)
    w = b*c
    e = fma(-b, c, w)
    f = fma(a, d, -w)
    return f+e
end


function area_lo2hi(a,b,c)
    p = (c + (b + a))*(a - (c - b))*(a + (c - b))*(c + (b - a))
    return sqrt(p)/4
end

function tantheta(pt1::NTuple{3,T}, pt2::NTuple{3,T}) where {T<:Real}
    h1, h2 = minmax(hypot(pt1...), hypot(pt2...))
    h3 = hypot((pt2 .- pt1)...)
    a  = 2*area_lo2hi(h1, h2, h3)
    return a/dot(p1,p2)
end
        
"""
    angle(point1::T, point2::T) where {T}

Accurately ascertains the undirected angle (0 <= radians < pi)
between two points given in N-dimensional Cartesian coordinates.

Prefer this to `acos` alternatives
- more reliably accurate
- more consistently stable

Suggested when any |coordinate| of either point may be outside 2^±20 or [1/1_000_000, 1_000_000].
Strongly recommended when any |coordinate| is outside 2^±24 or [1/16_000_000, 16_000_000].

If one of the points is at the origin, the result is zero.

You *must* define a tuple constructor `Tuple(x::YourPointType) = ...` if one does not already exist.
"""
function angle(point1::A, point2::A) where {N,T<:Real,NT<:NTuple{N,T}, V<:Vector{T}, A<:Union{NT,V}}
    unitpoint1 = unitize(point1)
    unitpoint2 = unitize(point2)

    y = unitpoint1 .- unitpoint2
    x = unitpoint1 .+ unitpoint2

    a = 2 * atan(norm(y) / norm(x))

    !(signbit(a) || signbit(T(pi) - a)) ? a : (signbit(a) ? zero(T) : T(pi))
end

@inline angle(point1::T, point2::T) where {T} = angle(Tuple(point1), Tuple(point2))

end # AngleBetweenVectors
