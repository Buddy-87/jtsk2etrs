import .Base: sin, cos, tan, asin,
             acos, atan, sqrt, max, min, minmax, ^

include("./ellipsoid.jl")

import .ell
"""
Converting cartesian perpendicular cartesian coordinates x,y,z
to geodetic latitude, longtitude, ellipsoidal height.
# Arguments
- `x::T` x-coordinate in [m]
- `y::T` y-coordinate in [m]
- `z::T` z-coordinate in [m]
- `ell::Ellipsoid` reference ellipsoid

# Returns
- `b::T` geodetic latitude in degrees
- `l::T` geodetic longitude in degrees
- `h::T` ellipsoid height in [m]

# Examples
```jldoctest
julia> xyz2blh(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6, wgs84)
(49.20589174, 16.59283443, 324.272) 
```

"""
function xyz2blh(x::T, y::T, z::T, ell::ellipsoid.Ellipsoid) where (T<:Real)
    l::T = atan(y,x)

    if l < 0 
        l+= 2.0*pi
    end

    p::T = sqrt(x^2 + y^2)
    b::T = atan(z/(p*(1.0 - ell.eccentricity)))
    nn::T = ell.a_axis/sqrt(1.0 - ell.eccentricity * (sin(b)^2));
    h::T = p/cos(b) - nn;

    db::T = 999
    dh::T = 999

    counter::Int8 = 0

    while  (db>1e-12 && dh>1e-6 || counter != 10)
        bi::T = atan((z*(nn+h))/(p*(h+nn*(1.0-ell.eccentricity))))
        Ni::T = ell.a_axis/sqrt(1.0 - ell.eccentricity * (sin(b))^2)
        hi::T = p/cos(b) - Ni

        db = abs(b-bi); dh = abs(h-hi);
        b = bi; h = hi; nn = Ni;
        counter+=1
    end

    return rad2deg(b), rad2deg(l), h
end

"""
Converting geodetic ellipsoidal coordinates b,l,h
to cartesian x, y, z.
# Arguments
- `b::T` geodetic latitude in degrees
- `l::T` geodetic longitude in degrees
- `h::T` ellipsoid height in [m]
- `ell::Ellipsoid` reference ellipsoid

# Returns
- `x::T` x-coordinate in [m]
- `y::T` y-coordinate in [m]
- `z::T` z-coordinate in [m]

# Examples
```jldoctest
julia> blh2xyz( 49.20589174, 16.59283443, 324.272, wgs84)
(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6) 
```

"""
function blh2xyz(b::T, l::T, h::T, ell::ellipsoid.Ellipsoid) where (T<:Real)
    brad::T = deg2rad(b)
    lrad::T = deg2rad(l)

    n::T = ell.a_axis/sqrt(1.0 - ell.eccentricity * (sin(brad))^2)
    x::T = (n+h)*cos(brad)*cos(lrad)
    y::T = (n+h)*cos(brad)*sin(lrad)
    z::T = (n*(1.0-ell.eccentricity) +h) * sin(brad)

    return x,y,z
end

"""
Converting the geodetic coordinates b,l,h to cartesian perpendicular cartesian coordinates x,y,z.
# Arguments
- `b::Vector{T}` column vector of geodetic latitude in degrees
- `l::Vector{T}` column vector of geodetic longitude in degrees
- `h::Vector{T}` column vector of ellipsoid height in [m]
- `ell::Ellipsoid` reference ellipsoid
# Returns
- `x::Vector{T}` column vector of x-coordinate in [m]
- `y::Vector{T}` column vector of y-coordinate in [m]
- `z::Vector{T}` column vector of z-coordinate in [m]
"""
function blh2xyz(b::Vector{T}, l::Vector{T}, h::Vector{T}, ell::ellipsoid.Ellipsoid) where (T<:Real)
    brad::Vector{T} = deg2rad.(b)
    lrad::Vector{T} = deg2rad.(l)

    n::Vector{T} = ell.a_axis/sqrt.(1.0 - ell.eccentricity .* (sin.(brad)).^2)
    x::Vector{T} = (n+h).*cos.(brad).*cos.(lrad)
    y::Vector{T} = (n+h).*cos.(brad).*sin.(lrad)
    z::Vector{T} = (n.*(1.0-ell.eccentricity) +h) .* sin.(brad)

    return x,y,z
end

#= function xyz2blh(x::Vector{T}, y::Vector{T}, z::Vector{T}, ell::ellipsoid.Ellipsoid) where (T<:Real)
    l::Vector{T} = atan.(y,x)

    if l < 0 
        l+= 2.0*pi
    end

    p::Vector{T} = sqrt.(x.^2 + y.^2)
    b::Vector{T} = atan.(z/(p.*(1.0 - ell.eccentricity)))
    nn::Vector{T} = ell.a_axis/sqrt.(1.0 - ell.eccentricity .* (sin.(b)).^2)
    h::Vector{T} = p/cos.(b) - nn

    db::Vector{T} = 999
    dh::Vector{T} = 999

    counter::Int8 = 0

    while  =#