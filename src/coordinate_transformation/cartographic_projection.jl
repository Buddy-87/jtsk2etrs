import .Base: sin, cos, tan, asin, deg2rad, rad2deg, abs
             acos, atan, sqrt, max, min, minmax, ^

"""
Converting the geodetic latitude and longitude (reference ellipsoid is Bessel)
to 2D cartographic coordinates in the JTSK projection. The transformation takes
multiple steps (workflow):
- geodetic coordinates on Bessel's ellipsoid (b,lambda)
- spherical coordinates on reference sphere (U,V)
- cartographic coordinates on reference sphere with shifted pole (S,D)
- polar coordinates on the cone (rho, epsilon) 
- rectangular coordinates in the plane (JTSK coordinates) (y, x)

Coordinates of the shifted pole are phi_Q = 48°15', lambda_Q = 42°30'.
# Arguments
- `b::T` geodetic latitude in DEG on Bessel's ellipsoid
- `l::T` geodetic longitude in DEG on Bessel's ellipsoid

# Returns 
- `y::T` y-coordinate in JTSK reference frame, +y-axis is in the west direction
- `x::T` x-coordinate in JTSK reference frame, +x-axis is in the south direction


# Example
```jldoctest
julia> bl2jtsk(,)
(,) 
```
"""
function bl2jtsk(b::T, l::T) where (T<:Real)
    const rho_0 = 1298039.004638987
    const s_0   = deg2rad(78.50000)
    const n     = 0.9799247046208300
    const alpha = 1.000597498371542
    const k     = 0.996592486900000 # 1.003419163966575
    const e     = 0.816968310215303e-1
    const vk    = deg2rad(42. + 31.0/60. + 31.41725/3600.0)
    const uk    = deg2rad(59. + 42.0/60. + 42.6968885/3600.0)
    const u0    = deg2rad(49. + 27.0/60. + 32.84625/3600.0)

    # (b,l) >> (U,V)
    phi::T = deg2rad(b);
    u::T = 0.2e1 * atan( pow(tan(phi/ 2.0 + pi/4.0 )
                                 * pow( (1. - e * sin(phi))
                                      / (1. + e * sin(phi)), e / 2.), alpha) / k)
               - pi/2.0;
    v::T = alpha * deg2rad(l + 17+36.0/60.0+46.002/3600.0); # CUZK -> 17°40'

    # (U,V) >> (S,D)
    s::T = asin(sin(u) * sin(uk) + cos(u) * cos(uk) * cos(vk - v));
    d::T = asin(cos(u) * sin(vk - v) / cos(s));

    # (S,D) >> (rho, epsilon)
    rho::T = rho_0 * (tan(s_0 /2. + pi/4.0) / tan(s / 2.0 + pi/4.0))^n
    epsilon::T = n * d;

    # (\rho, \vaerpsilon) >> (Y,X)
    x::T = rho * cos(epsilon);
    y::T = rho * sin(epsilon);

    return y, x;
end

"""
Converting 2D cartographic coordinates from the JTSK projection (y,x) 
to the geodetic latitude and longitude (reference ellipsoid is Bessel). 
The transformation takes multiple steps (workflow):
- rectangular coordinates in the plane (JTSK coordinates) (y, x)
- polar coordinates on the cone (rho, epsilon) 
- cartographic coordinates on reference sphere with shifted pole (S,D)
- spherical coordinates on reference sphere (U,V)
- geodetic coordinates on Bessel's ellipsoid (b,lambda)

Coordinates of the shifted pole are phi_Q = 48°15', lambda_Q = 42°30'.
# Arguments
- `y::T` y-coordinate in JTSK reference frame, +y-axis is in the west direction
- `x::T` x-coordinate in JTSK reference frame, +x-axis is in the south direction

# Returns 
- `b::T` geodetic latitude in DEG on Bessel's ellipsoid
- `l::T` geodetic longitude in DEG on Bessel's ellipsoid


# Example
```jldoctest
julia> jtsk2bl(,)
(,) 
```
"""
function jtsk2bl( y::T, x::T) where (T<:Real)
    const rho_0 = 1298039.004638987;
    const s_0   = deg2rad(78.50000) ;
    const n     = .9799247046208300;
    const alpha = 1.000597498371542;
    const k     = /*1.003419163966575;*/.996592486900000;
    const e     = 0.81696831215303e-1;
    const vk    = deg2rad(42. + 31.0/60. + 31.41725/3600.0)
    const uk    = deg2rad(59. + 42.0/60. + 42.6968885/3600.0)
    const u0    = deg2rad(49. + 27.0/60. + 32.84625/3600.0)

    # ((\rho, \vaerpsilon) << (Y,X)
    rho::T = sqrt(x*x+y*y);
    epsilon::T = atan(y/x);

    # (S,D) << (\rho, \vaerpsilon)
    d::T = epsilon/n ;
    s::T = 2.0*atan((rho_0/rho)^(1.0/n) * tan( s_0/2. + pi/4.0)) - pi/2.0;

    # (U,V) << (S,D)
    u::T = asin( sin(uk)*sin(s) - cos(uk)*cos(s)*cos(d));
    delta_v::T = asin( sin(d)*cos(s)/cos(u));

    # (phi,lambda) << (U,V)
    lambda::T = (vk - delta_v)/alpha - 0.3074009723405860142863407 ;
    phi_i::T = u 
    phi_ii::T = u
    dphi::T =9999
    loop_cnt::Int8 = 0

    while ( abs( dphi ) >= .000001/3600*DEG2RAD )
        phi_ii = 2.0* atan((k * tan( u/2.0 + pi/4.0 ))^(1.0/alpha) 
                         * ((1.0-e*sin(phi_i))/(1.0+e*sin(phi_i)))^(-e/2.0 )) - pi/2.0;
        dphi = phi_ii - phi_i;
        phi_i = phi_ii;

        if (loop_cnt > 10 )
            break
        end

        loop_cnt += 1; # approx 5 iterations
    end

    b::T = rad2deg(phi_ii)
    l::T = rad2deg(lambda)

    return b,l
end


"""
Converting the geodetic latitude and longitude (reference ellipsoid is Bessel)
to 2D cartographic coordinates in the JTSK05 projection. The transformation takes
multiple steps (workflow):
- geodetic coordinates on Bessel's ellipsoid (b,lambda)
- spherical coordinates on reference sphere (U,V)
- cartographic coordinates on reference sphere with shifted pole (S,D)
- polar coordinates on the cone (rho, epsilon) 
- rectangular coordinates in the plane (JTSK coordinates) (y, x)
- add corrections from cubic transformation

Coordinates of the shifted pole are phi_Q = 48°15', lambda_Q = 42°30'.
# Arguments
- `b::T` geodetic latitude in DEG on Bessel's ellipsoid
- `l::T` geodetic longitude in DEG on Bessel's ellipsoid

# Returns 
- `y::T` y-coordinate in JTSK reference frame, +y-axis is in the west direction
- `x::T` x-coordinate in JTSK reference frame, +x-axis is in the south direction


# Example
```jldoctest
julia> bl2jtsk(,)
(,) 
```
"""
function bl2jtsk05( b::T, l::T ) where (T<:Real)
    y::T,x::T = bl2jtsk::T(b, l)
    dy::T, dx::T = interpolate_correction(y,x);

    return y - dy, x - dx;
end

"""
Converting 2D cartographic coordinates from the JTSK05 projection (y,x) 
to the geodetic latitude and longitude (reference ellipsoid is Bessel). 
The transformation takes multiple steps (workflow):
- substract corrections from cubic transformation
- rectangular coordinates in the plane (JTSK coordinates) (y, x)
- polar coordinates on the cone (rho, epsilon) 
- cartographic coordinates on reference sphere with shifted pole (S,D)
- spherical coordinates on reference sphere (U,V)
- geodetic coordinates on Bessel's ellipsoid (b,lambda)

Coordinates of the shifted pole are phi_Q = 48°15', lambda_Q = 42°30'.
# Arguments
- `y::T` y-coordinate in JTSK reference frame, +y-axis is in the west direction
- `x::T` x-coordinate in JTSK reference frame, +x-axis is in the south direction

# Returns 
- `b::T` geodetic latitude in DEG on Bessel's ellipsoid
- `l::T` geodetic longitude in DEG on Bessel's ellipsoid


# Example
```jldoctest
julia> jtsk052bl(,)
(,) 
```
"""
function jtsk052bl( y::T, x::T) where (T<:Real)
    dy::T, dx::T = interpolate_correction(y,x);
    return jtsk2bl(y + dy, x + dx);
end

"""

# Arguments
- `y::T` y-coordinate in JTSK reference frame, +y-axis is in the west direction
- `x::T` x-coordinate in JTSK reference frame, +x-axis is in the south direction

# Returns 
- `dy::T` correction for transforming the JTSK to JTSK05 in y-axis
- `dx::T` correction for transforming the JTSK to JTSK05 in x-axis

# Example
```jldoctest
julia> interpolate_correction(,)
(,) 
```
"""
function interpolate_correction(y::T, x::T) where (T<:Real)
    # Cubic transformation constants
    const a1  =  0.2946529277e-01
    const a2  =  0.2515965696e-01
    const a3  =  0.1193845912e-06
    const a4  = -0.4668270147e-06
    const a5  =  0.9233980362e-11
    const a6  =  0.1523735715e-11
    const a7  =  0.1696780024e-17
    const a8  =  0.4408314235e-17
    const a9  = -0.8331083518e-23
    const a10 = -0.3689471323e-23

    x_red::T = x - 1089000
    y_red::T = y - 654000

    dy::T = a2+a3*y_red + a4*x_red + 2.0*a5*y_red*x_red + a6*(x_red*x_red - y_red*y_red)
            +a8*x_red*(x_red*x_red - 3.0*y_red*y_red) + a7*y_red(3.0*x_red*x_red - y_red*y_red)
            -4.0*a10*x_red*y_red*(x_red*x_red - y_red*y_red)
            +a9*( x_red^4.0 + y_red^4.0 - 6.0*(x_red*y_red)^2.0)

    dx::T = a1+a3*x_red-a4*y_red - 2.0*a6*x_red*y_red + a5*( x_red*x_red - y_red*y_red )
            +a7*x_red*(x_red*x_red-3.0*y_red*y_red) - a8*y_red(3.0*x_red*x_red-y_red*y_red )
            +4.0*a9*y_red*x_red(x_red*x_red-y_red*y_red)
            +10.0*( x_red^4.0 +y_red^4.0 - 6.0*(x_red*y_red)^2.0 )

    return dy, dx
end