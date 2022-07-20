import .Base: sin, cos, tan, asin,
             acos, atan, sqrt, max, min, minmax, ^,
             cosh, sinh, asinh, acosh, rad2deg, deg2rad


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
function xyz2blh(x::T, y::T, z::T, ell::Ellipsoid) where (T<:AbstractFloat)
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
function blh2xyz(b::T, l::T, h::T, ell::Ellipsoid) where (T<:AbstractFloat)
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
function blh2xyz(b::Vector{T}, l::Vector{T}, h::Vector{T}, ell::Ellipsoid) where (T<:AbstractFloat)
    brad::Vector{T} = deg2rad.(b)
    lrad::Vector{T} = deg2rad.(l)

    n::Vector{T} = ell.a_axis./sqrt.(1.0 .- ell.eccentricity .* (sin.(brad)).^2)
    x::Vector{T} = (n+h).*cos.(brad).*cos.(lrad)
    y::Vector{T} = (n+h).*cos.(brad).*sin.(lrad)
    z::Vector{T} = (n.*(1.0-ell.eccentricity) .+ h) .* sin.(brad)

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
function xyz2blh(x::Vector{T}, y::Vector{T}, z::Vector{T}, ell::Ellipsoid) where (T<:AbstractFloat)
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
    b::Vector{T} = atan.(z./(p.*(1.0 .- ell.eccentricity)))
    nn::Vector{T} = ell.a_axis./sqrt.(1.0 .- ell.eccentricity .* (sin.(b)).^2)
    h::Vector{T} = p./cos.(b) .- nn

    counter::Int8 = 0

    while true
        bi::Vector{T} = atan.((z.*(nn+h))./(p.*(h+nn.*(1.0.-ell.eccentricity))))
        Ni::Vector{T} = ell.a_axis./sqrt.(1.0 .- ell.eccentricity .* (sin.(b)).^2)
        hi::Vector{T} = p./cos.(b) .- Ni 
        # hi::T = p/cos(b) - Ni

        b = bi; h = hi; nn = Ni;
        counter+=1

        if ( maximum(abs.(bi .- b)) > 1e-10  || counter == 12)
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
    phi::T = atan(y,x)

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
    phi::Vector{T} = atan.(y,x)

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
end

"""
Converts the cartesian coordinates x,y,z to ellipsoidal curvilinear coordinates u,v,w.
Where the relation between the cartesian coordinates x,y,z and the ellipsoidal coordinates u,v,w is
- `x = a e sin(u) cos(v) cosh(w)`
- `y = a e sin(u) sin(v) cosh(w)`
- `z = a e cos(u) sinh(w)`
# Arguments
- `x::T` x-coordinate in [m]
- `y::T` y-coordinate in [m]
- `z::T` z-coordinate in [m]
# Returns
- `u::T` u-coordinate in DEG format
- `v::T` v-coordinate in DEG format
- `w::T` w-coordinate is unit-less
# Examples
```jldoctest
julia> xyz2uvw(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6, wgs84)
(44.972761953232066, 15.987500000000004, 3.1947638833052543)
```
"""
function xyz2uvw(x::T, y::T, z::T, ell::Ellipsoid) where (T<:AbstractFloat)
    a::T = ell.a_axis
    e::T = sqrt(ell.eccentricity)
    ae::T = a*e
    v::T = atan(y,x)

    p::T = sqrt(x^2 + y^2)
    r::T = sqrt(p^2 + z^2)

    s::T = sqrt((ae^2.0 + r*r)^2.0 - 4.0*ae^2.0*p*p)

    dd1::T = ((ae)^(2.0) + r*r + s)/(2.0*(ae)^(2.0))
    dd2::T = ((ae)^(2.0) + r*r - s)/(2.0*(ae)^(2.0))

    dd::T = 1.0
    if ( abs(dd1) > 1.0 && dd2 == 1.0) 
        dd = dd2 - 1.0e-10
    elseif ( abs(dd1) > 1.0 && dd1 != 1.0)
        dd = dd2
    else
        dd = dd1
    end

    f::T = z/(ae*sqrt(1.0 - dd))
    w::T = 1.0
    u::T = 1.0
    if z > 0.0
        w = asinh(f)
        u = asin(sqrt(dd))
    elseif z < 0.0
        w = asinh(-f)
        u = pi - asin(sqrt(dd)) 
    else
        w = acosh(1.0/e)
        u = asin(sqrt(dd-1.0e-10));
    end

    return rad2deg(u), rad2deg(v), w
end

"""
Converts the cartesian coordinates x,y,z to ellipsoidal curvilinear coordinates u,v,w.
Where the relation between the cartesian coordinates x,y,z and the ellipsoidal coordinates u,v,w is
- `x = a e sin(u) cos(v) cosh(w)`
- `y = a e sin(u) sin(v) cosh(w)`
- `z = a e cos(u) sinh(w)`
# Arguments
- `x::Vector{T}` vector with x-coordinate in [m]
- `y::Vector{T}` vector with y-coordinate in [m]
- `z::Vector{T}` vector with z-coordinate in [m]
# Returns
- `u::Vector{T}` vector with u-coordinate in DEG format
- `v::Vector{T}` vector with v-coordinate in DEG format
- `w::Vector{T}` vector with w-coordinate is unit-less

"""
function xyz2uvw(x::Vector{T}, y::Vector{T}, z::Vector{T}, ell::Ellipsoid) where (T<:AbstractFloat)
    a::T = ell.a_axis
    e::T = sqrt(ell.eccentricity)
    ae::T = a*e
    v::Vector{T} = atan.(y,x)

    p::Vector{T} = sqrt.(x.^2 + y.^2)
    r::Vector{T} = sqrt.(p.^2 + z.^2)

    s::Vector{T} = sqrt.((ae^2.0 + r.*r)^2.0 - 4.0*ae^2.0*p.*p)

    dd1::Vector{T} = ((ae)^(2.0) + r.*r + s)./(2.0*(ae)^(2.0))
    dd2::Vector{T} = ((ae)^(2.0) + r.*r - s)./(2.0*(ae)^(2.0))

    dd::Vector{T} = ones(size(x))

    for i=1:length(x)
        if ( abs(dd1[i]) > 1.0 && dd2[i] == 1.0) 
            dd[i] = dd2[i] - 1.0e-10
        elseif ( abs(dd1[i]) > 1.0 && dd1[i] != 1.0)
            dd[i] = dd2[i]
        else
            dd[i] = dd1[i]
        end
    end

    f::Vector{T} = z./(ae*sqrt.(1.0 - dd))
    w::Vector{T} = ones(size(x))
    u::Vector{T} = ones(size(x))

    for i=1:length(x)
        if z[i] > 0.0
            w[i] = asinh(f[i])
            u[i] = asin(sqrt(dd[i]))
        elseif z[i] < 0.0
            w[i] = asinh(-f[i])
            u[i] = pi - asin(sqrt(dd[i])) 
        else
            w[i] = acosh(1.0/e)
            u[i] = asin(sqrt(dd[i]-1.0e-10));
        end
    end

    return rad2deg.(u), rad2deg.(v), w
end

"""
Converts the cartesian coordinates x,y,z to ellipsoidal curvilinear coordinates u,v,w.
Where the relation between the cartesian coordinates x,y,z and the ellipsoidal coordinates u,v,w is
- `x = a e sin(u) cos(v) cosh(w)`
- `y = a e sin(u) sin(v) cosh(w)`
- `z = a e cos(u) sinh(w)`
# Arguments
- `u::T` u-coordinate in DEG format
- `v::T` v-coordinate in DEG format
- `w::T` w-coordinate is unit-less
# Returns
- `x::T` x-coordinate in [m]
- `y::T` y-coordinate in [m]
- `z::T` z-coordinate in [m]
# Examples
```jldoctest
uvw2xyz(44.972761953232066, 15.987500000000004, 3.1947638833052543, wgs84)
(4.333743281298126e6, 1.241657737047367e6, 4.49726939942041e6)
```
"""
function uvw2xyz(u::T, v::T, w::T, ell::Ellipsoid) where (T<:AbstractFloat)
    a::T = ell.a_axis
    e::T = sqrt(ell.eccentricity)
    ae::T = a*e

    uu::T = deg2rad(u)
    vv::T = deg2rad(v)

    x::T = ae*sin(uu)*cos(vv)*cosh(w)
    y::T = ae*sin(uu)*sin(vv)*cosh(w)
    z::T = ae*cos(uu)*sinh(w)

    return x,y,z
end

    
"""
Converts the cartesian coordinates x,y,z to ellipsoidal curvilinear coordinates u,v,w.
Where the relation between the cartesian coordinates x,y,z and the ellipsoidal coordinates u,v,w is
- `x = a e sin(u) cos(v) cosh(w)`
- `y = a e sin(u) sin(v) cosh(w)`
- `z = a e cos(u) sinh(w)`
# Arguments
- `u::Vector{T}` u-coordinate in DEG format
- `v::Vector{T}` v-coordinate in DEG format
- `w::Vector{T}` w-coordinate is unit-less
# Returns
- `x::Vector{T}` x-coordinate in [m]
- `y::Vector{T}` y-coordinate in [m]
- `z::Vector{T}` z-coordinate in [m]
# Examples
```jldoctest
julia> uvw2xyz(44.972761953232066, 15.987500000000004, 3.1947638833052543, wgs84)
(4.333743281298126e6, 1.241657737047367e6, 4.49726939942041e6)
```
"""
function uvw2xyz(u::Vector{T}, v::Vector{T}, w::Vector{T}, ell::Ellipsoid) where (T<:AbstractFloat)
    a::T = ell.a_axis
    e::T = sqrt(ell.eccentricity)
    ae::T = a*e

    uu::Vector{T} = deg2rad.(u)
    vv::Vector{T} = deg2rad.(v)

    x::Vector{T} = ae*sin(uu).*cos(vv).*cosh(w)
    y::Vector{T} = ae*sin(uu).*sin(vv).*cosh(w)
    z::Vector{T} = ae*cos(uu).*sinh(w)

    return x,y,z

end