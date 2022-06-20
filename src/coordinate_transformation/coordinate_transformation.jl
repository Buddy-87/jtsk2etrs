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
function xyz2blh(x::T, y::T, z::T, ell::ellipsoid.Ellipsoid) where (T<:AbstractFloat)
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
function blh2xyz(b::T, l::T, h::T, ell::ellipsoid.Ellipsoid) where (T<:AbstractFloat)
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
function blh2xyz(b::Vector{T}, l::Vector{T}, h::Vector{T}, ell::ellipsoid.Ellipsoid) where (T<:AbstractFloat)
    brad::Vector{T} = deg2rad.(b)
    lrad::Vector{T} = deg2rad.(l)

    n::Vector{T} = ell.a_axis/sqrt.(1.0 - ell.eccentricity .* (sin.(brad)).^2)
    x::Vector{T} = (n+h).*cos.(brad).*cos.(lrad)
    y::Vector{T} = (n+h).*cos.(brad).*sin.(lrad)
    z::Vector{T} = (n.*(1.0-ell.eccentricity) +h) .* sin.(brad)

    return x,y,z
end

"""
Converting cartesian perpendicular cartesian coordinates x,y,z
to geodetic latitude, longtitude, ellipsoidal height.
# Arguments
- `x::Vector{T}` column vector of x-coordinate in [m]
- `y::Vector{T}` column vector of y-coordinate in [m]
- `z::Vector{T}` column vector of z-coordinate in [m]
- `ell::Ellipsoid` reference ellipsoid
# Returns
- `b::Vector{T}` column vector of geodetic latitude in degrees
- `l::Vector{T}` column vector of geodetic longitude in degrees
- `h::Vector{T}` column vector of ellipsoid height in [m]
"""
function xyz2blh(x::Vector{T}, y::Vector{T}, z::Vector{T}, ell::ellipsoid.Ellipsoid) where (T<:AbstractFloat)
    l::Vector{T} = atan.(y,x)

    # Find all negative values of l and add 2pi to them
    # zindex = findall(<(0.0), l)
    # l[CartesianIndex.(zindex)] .+= 2.0*pi
    for i in 1:length(l)
        if l[i] < 0
            l[i] += 2.0*pi
        end
    end

    p::Vector{T} = sqrt.(x.^2 + y.^2)
    b::Vector{T} = atan.(z/(p.*(1.0 - ell.eccentricity)))
    nn::Vector{T} = ell.a_axis/sqrt.(1.0 - ell.eccentricity .* (sin.(b)).^2)
    h::Vector{T} = p/cos.(b) - nn

    db::Vector{T} = 999
    dh::Vector{T} = 999

    counter::Int8 = 0

    while true
        bi::Vector{T} = atan.((z.*(nn+h))./(p.*(h+nn.*(1.0-ell.eccentricity))))
        Ni::Vector{T} = ell.a_axis/sqrt.(1.0 - ell.eccentricity .* (sin.(b)).^2)
        hi::Vector{T} = p/cos.(b) - Ni

        b = bi; h = hi; nn = Ni;
        counter+=1

        if ( maximum(abs.( bi .- b))  || counter == 10)
            break
        end
    end

    return rad2deg.(b), rad2deg.(l), h
end

"""
Converting the cartesian coordinates x,y,z to spherical coordinates r,theta,phi.
# Arguments
- `x::T` x-coordinate in [m]
- `y::T` y-coordinate in [m]
- `z::T` z-coordinate in [m]
# Returns
- `r::T` radius in [m]
- `theta::T` polar angle in degrees
- `phi::T` azimuthal angle in degrees
# Examples
```jldoctest
julia> xyz2sph(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6)
(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6)
```
"""
function xyz2sph(x::T, y::T, z::T) where (T<:AbstractFloat)
    r::T = sqrt(x^2 + y^2 + z^2)
    theta::T = atan(z/sqrt(x^2 + y^2))
    phi::T = atan(y/x)

    if phi < 0 
        phi+= 2.0*pi
    end

    return r, rad2deg(theta), rad2deg(phi)
end

"""
Converting the cartesian coordinates x,y,z to spherical coordinates r,theta,phi.
# Arguments
- `x::Vector{T}` column vector of x-coordinate in [m]
- `y::Vector{T}` column vector of y-coordinate in [m]
- `z::Vector{T}` column vector of z-coordinate in [m]
# Returns
- `r::Vector{T}` column vector of radius in [m]
- `theta::Vector{T}` column vector of polar angle in degrees
- `phi::Vector{T}` column vector of azimuthal angle in degrees
# Examples
```jldoctest
julia> xyz2sph(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6)
(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6)
```
"""
function xyz2sph(x::Vector{T}, y::Vector{T}, z::Vector{T}) where (T<:AbstractFloat)
    r::Vector{T} = sqrt.(x.^2 + y.^2 + z.^2)
    theta::Vector{T} = atan.(z./sqrt.(x.^2 + y.^2))
    phi::Vector{T} = atan.(y./x)

    for i in 1:length(phi)
        if phi[i] < 0
            phi[i] += 2.0*pi
        end
    end

    return r, rad2deg.(theta), rad2deg.(phi)
end

"""
Converting the cartesian coordinates x,y,z to spherical coordinates r,theta,phi.
# Arguments
- `r::T` radius in [m]
- `theta::T` polar angle in degrees
- `phi::T` azimuthal angle in degrees
# Returns
- `x::T` x-coordinate in [m]
- `y::T` y-coordinate in [m]
- `z::T` z-coordinate in [m]
# Examples
```jldoctest
julia> sph2xyz(12.16,47.416,112.231)
(-3.1131133109015656, 7.616668705952813, 8.953238706593643)
```
"""
function sph2xyz(r::T, theta::T, phi::T) where (T<:AbstractFloat)
    u::T = deg2rad(theta)
    v::T = deg2rad(phi)

    x::T = r*cos(u)*cos(v)
    y::T = r*cos(u)*sin(v)
    z::T = r*sin(u)

    return x,y,z
end

"""
Converting the cartesian coordinates x,y,z to spherical coordinates r,theta,phi.
# Arguments
- `x::Vector{T}` column vector of x-coordinate in [m]
- `y::Vector{T}` column vector of y-coordinate in [m]
- `z::Vector{T}` column vector of z-coordinate in [m]
# Returns
- `r::Vector{T}` column vector of radius in [m]
- `theta::Vector{T}` column vector of polar angle in degrees
- `phi::Vector{T}` column vector of azimuthal angle in degrees
# Examples
```jldoctest
julia> xyz2sph(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6)
(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6)
```
"""
function sph2xyz(r::Vector{T}, theta::Vector{T}, phi::Vector{T}) where (T<:AbstractFloat)
    u::Vector{T} = deg2rad.(theta)
    v::Vector{T} = deg2rad.(phi)

    x::Vector{T} = r.*cos.(u).*cos.(v)
    y::Vector{T} = r.*cos.(u).*sin.(v)
    z::Vector{T} = r.*sin.(u)

    return x,y,z
end;