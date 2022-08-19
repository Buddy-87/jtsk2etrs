import .Base: sin, cos, tan, asin, deg2rad, rad2deg, abs
             acos, atan, sqrt, max, min, minmax, ^
"""
TUBO 49° 12' 21.20976"	 	16° 35' 34.20421"	 	324.374	 	0.3107
     49.2058916             16.592834502777777
     599 131.597            1 159 442.044

bTubo = 49.2058916;
lTubo = 16.592834502777777;
hTubo = 324.374;

xt, yt, zt = blh2xyz(bTubo, lTubo, hTubo, wgs84)
xbt,ybt,zbt=etrf2bessel(xt,yt,zt)
bbt, lbt, hbt = xyz2blh(xbt, ybt, zbt, bessel)

(49.204101717107235, 16.591509006826552, 1497.4366400362924)

bl2jtsk(49.204101717107235, 16.591509006826552)
"""


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

Constants :

``\\varphi_0 = 49^\\circ 30'``

``\\alpha = \\sqrt{1 + \\frac{e^2 \\cos^4 \\varphi_0}{1-e^2}}``

`` U_0 = \\asin \\frac{\\sin \\varphi_0}{\\alpha} ``

`` g_{\\varphi_0} = \\left( \\frac{1 + e \\sin \\varphi_0}{1 - e \\sin \\varphi_0} \\right)^{\\alpha e / 2} ``

# Arguments
- `b::T` geodetic latitude in DEG on Bessel's ellipsoid
- `l::T` geodetic longitude in DEG on Bessel's ellipsoid

# Returns 
- `y::T` y-coordinate in JTSK reference frame, +y-axis is in the west direction
- `x::T` x-coordinate in JTSK reference frame, +x-axis is in the south direction


# Example
```jldoctest
# b = 50˚ 57' 8.3936'', l = 14° 34' 51.1547''
julia> bl2jtsk(50.95233155555556, 14.580876305555556)
(,) 
```
"""
function bl2jtsk(b::T, l::T) where (T<:AbstractFloat)
    rho_0::T = 1298039.00463898700
    s_0::T   = 1.37008346281554870 # 78°30'
    n0::T    = 0.97992470462083000
    alpha::T = 1.00059749837154240 # s
    k::T     = 1.00341916396657260 # 0.99659248690000000
    e::T     = 0.08169683121529256
    g::T     = 1.00509776241540980 # `` g( \varphi_0) ``
    vk::T    = 0.74220813543248410 # 42°31'31.41725'' 
    uk::T    = 1.04216856380474300 # 59˚42'42.6968885''
    u0::T    = 0.86323910265848790

    # (b,l) >> (U,V)
    phi::T = deg2rad(b);
    gphi::T = ((1.0 + e * sin(phi))/(1.0 - e * sin(phi)))^(alpha * e / 2.0)

    # u::T = 2.0 * (atan(k * (tan(1.2173671532660448))^(alpha)/gphi) - pi/4.0)
    u::T = 2.0 * (atan(k * (tan(phi/2.0 + pi/4.0))^(alpha)/gphi) - pi/4.0)
    dv::T = alpha * deg2rad(24.833333333333332 - l)

    # (U,V) >> (S,D)
    s::T = asin(sin(u) * sin(uk) + cos(u) * cos(uk) * cos(dv));
    d::T = asin(cos(u) * sin(dv) / cos(s));

    # (S,D) >> (rho, epsilon)
    rho::T = rho_0 * (tan(s_0 /2.0 + pi/4.0) / tan(s / 2.0 + pi/4.0))^n0
    epsilon::T = n0 * d;

    # (\rho, \vaerpsilon) >> (Y,X)
    x::T = rho * cos(epsilon);
    y::T = rho * sin(epsilon);

    return y, x;
end

function bl2jtsk(b::Vector{T}, l::Vector{T}) where (T<:AbstractFloat)
    rho_0::T = 1298039.00463898700
    s_0::T   = 1.37008346281554870 # 78°30'
    n0::T    = 0.97992470462083000
    alpha::T = 1.00059749837154240 # s
    k::T     = 1.00341916396657260 # 0.99659248690000000
    e::T     = 0.08169683121529256
    g::T     = 1.00509776241540980 # `` g( \varphi_0) ``
    vk::T    = 0.74220813543248410 # 42°31'31.41725'' 
    uk::T    = 1.04216856380474300 # 59˚42'42.6968885''
    u0::T    = 0.86323910265848790

    # (b,l) >> (U,V)
    phi::Vector{T} = deg2rad.(b);
    gphi::Vector{T} = ((1.0 .+ e .* sin.(phi))./(1.0 .- e .* sin.(phi))).^(alpha * e / 2.0)

    u::Vector{T} = 2.0 .* (atan.(k .* (tan.(phi/2.0 .+ pi/4.0)).^(alpha)./gphi) .- pi/4.0)
    dv::Vector{T} = alpha .* deg2rad.(24.833333333333332 .- l)

    # (U,V) >> (S,D)
    s::Vector{T} = asin.(sin.(u) .* sin.(uk) .+ cos.(u) .* cos.(uk) .* cos.(dv));
    d::Vector{T} = asin.(cos.(u) .* sin.(dv) ./ cos.(s));

    # (S,D) >> (rho, epsilon)
    rho::Vector{T} = rho_0 .* (tan.(s_0 ./2.0 .+ pi/4.0) ./ tan.(s ./ 2.0 .+ pi/4.0)).^n0
    epsilon::Vector{T} = n0 .* d;

    # (\rho, \vaerpsilon) >> (Y,X)
    x::Vector{T} = rho .* cos.(epsilon);
    y::Vector{T} = rho .* sin.(epsilon);

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
function jtsk2bl(y::T, x::T) where (T<:AbstractFloat)
    rho_0::T = 1298039.00463898700
    s_0::T   = 1.37008346281554870 # 78°30'
    n0::T    = 0.97992470462083000 # sin(s_0)
    alpha::T = 1.00059749837154240 # s
    k::T     = 1.00341916396657260 # 0.99659248690000000
    e::T     = 0.08169683121529256
    g::T     = 1.00509776241540980 # `` g( \varphi_0) ``
    vk::T    = 0.74220813543248410 # 42°31'31.41725'' 
    uk::T    = 1.04216856380474300 # 59˚42'42.6968885''
    u0::T    = 0.86323910265848790

    # ((\rho, \vaerpsilon) << (Y,X)
    rho::T = sqrt(x*x+y*y);
    epsilon::T = atan(y/x);

    # (S,D) << (\rho, \vaerpsilon)
    d::T = epsilon/n0;
    s::T = 2.0*(atan((rho_0/rho)^(1.0/n0) * tan(s_0/2. + pi/4.0)) - pi/4.0);

    # (U,V) << (S,D)
    u::T = asin( sin(uk)*sin(s) - cos(uk)*cos(s)*cos(d));
    dv::T = asin( sin(d)*cos(s)/cos(u));

    # (phi,lambda) << (U,V)
    lambda::T = (0.433423430912-dv/alpha)
    phi_i::T = u 
    phi_ii::T = u
    dphi::T =9999
    loop_cnt::Int8 = 0

    while ( abs( dphi ) >= deg2rad(0.000001/3600) )
        gprime_phi::T = ((1.0 + e * sin(phi_i))/(1.0 - e * sin(phi_i)))

        phi_ii = 2.0*(atan(k^(-1.0/alpha) * (tan(u/2.0 + pi/4.0))^(1.0/alpha) * gprime_phi^(e/2.0) ) - pi/4.0)
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

function jtsk2bl(y::Vector{T}, x::Vector{T}) where (T<:AbstractFloat)
    n::Int64 = length(y)

    b::Vector{T} = Vector{T}(undef, n)
    l::Vector{T} = Vector{T}(undef, n)

    for i=1:n
        b[i], l[i] = jtsk2bl(y[i], x[i])
    end

    return b, l
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
function bl2jtsk05( b::T, l::T ) where (T<:AbstractFloat)
    y::T,x::T = bl2jtsk(b, l)
    dy::T, dx::T = interpolate_correction(y,x);

    return y - dy, x - dx;
end

function bl2jtsk05( b::Vector{T}, l::Vector{T} ) where (T<:AbstractFloat)
    y::Vector{T},x::Vector{T} = bl2jtsk(b, l)
    dy::Vector{T}, dx::Vector{T} = interpolate_correction(y,x);

    return y .- dy, x .- dx;
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
function jtsk052bl( y::T, x::T) where (T<:AbstractFloat)
    dy::T, dx::T = interpolate_correction(y,x);
    return jtsk2bl(y + dy, x + dx);
end

function jtsk052bl( y::Vector{T}, x::Vector{T}) where (T<:AbstractFloat)
    dy::Vector{T}, dx::Vector{T} = interpolate_correction(y,x);
    return jtsk2bl(y .+ dy, x .+ dx);
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
function interpolate_correction(y::T, x::T) where (T<:AbstractFloat)
    # Cubic transformation constants
    a1::T  =  0.2946529277e-01
    a2::T  =  0.2515965696e-01
    a3::T  =  0.1193845912e-06
    a4::T  = -0.4668270147e-06
    a5::T  =  0.9233980362e-11
    a6::T  =  0.1523735715e-11
    a7::T  =  0.1696780024e-17
    a8::T  =  0.4408314235e-17
    a9::T  = -0.8331083518e-23
    a10::T = -0.3689471323e-23

    x_red::T = x - 1089000.0
    y_red::T = y - 654000.0

    dy::T = a2+a3*y_red + a4*x_red + 2.0*a5*y_red*x_red + a6*(x_red*x_red - y_red*y_red)
            +a8*x_red*(x_red*x_red - 3.0*y_red*y_red) + a7*y_red*(3.0*x_red*x_red - y_red*y_red)
            -4.0*a10*x_red*y_red*(x_red*x_red - y_red*y_red)
            +a9*( x_red^4.0 + y_red^4.0 - 6.0*(x_red*y_red)^2.0)

    dx::T = a1+a3*x_red-a4*y_red - 2.0*a6*x_red*y_red + a5*( x_red*x_red - y_red*y_red )
            +a7*x_red*(x_red*x_red-3.0*y_red*y_red) - a8*y_red*(3.0*x_red*x_red-y_red*y_red)
            +4.0*a9*y_red*x_red*(x_red*x_red-y_red*y_red)
            +10.0*( x_red^4.0 +y_red^4.0 - 6.0*(x_red*y_red)^2.0 )

    return dy, dx
end

"""
Use cubic interpolation to compute the correction for transforming the JTSK to JTSK05.
These discrepences are caused by deformation of the original JTSK.
# Arguments
- `y::Vector{T}` y-coordinates in JTSK reference frame, +y-axis is in the west direction
- `x::Vector{T}` x-coordinates in JTSK reference frame, +x-axis is in the south direction
# Returns
- `dy::Vector{T}` correction for transforming the JTSK to JTSK05 in y-axis
- `dx::Vector{T}` correction for transforming the JTSK to JTSK05 in x-axis
# Example
```jldoctest
julia> interpolate_correction(,)
(,) 
```
""" 
function interpolate_correction(y::Vector{T}, x::Vector{T}) where (T<:AbstractFloat)
    # Cubic transformation constants
    a1::T  =  0.2946529277e-01
    a2::T  =  0.2515965696e-01
    a3::T  =  0.1193845912e-06
    a4::T  = -0.4668270147e-06
    a5::T  =  0.9233980362e-11
    a6::T  =  0.1523735715e-11
    a7::T  =  0.1696780024e-17
    a8::T  =  0.4408314235e-17
    a9::T  = -0.8331083518e-23
    a10::T = -0.3689471323e-23

    x_red::Vector{T} = x .- 1089000
    y_red::Vector{T} = y .- 654000

    dy::Vector{T} = a2.+a3.*y_red .+ a4.*x_red .+ 2.0.*a5.*y_red.*x_red .+ a6.*(x_red.*x_red .- y_red.*y_red)
            .+a8.*x_red.*(x_red.*x_red .- 3.0.*y_red.*y_red) .+ a7.*y_red.*(3.0.*x_red.*x_red .- y_red.*y_red)
            .-4.0.*a10.*x_red.*y_red.*(x_red.*x_red .- y_red.*y_red)
            .+a9.*( x_red.^4.0 .+ y_red.^4.0 .- 6.0.*(x_red.*y_red).^2.0)

    dx::Vector{T} = a1.+a3.*x_red-a4.*y_red .- 2.0.*a6.*x_red.*y_red .+ a5.*( x_red.*x_red .- y_red.*y_red )
            .+a7.*x_red.*(x_red.*x_red.-3.0.*y_red.*y_red) .- a8.*y_red.*(3.0.*x_red.*x_red.-y_red.*y_red)
            .+4.0.*a9.*y_red.*x_red.*(x_red.*x_red.-y_red.*y_red)
            .+10.0.*( x_red.^4.0 .+y_red.^4.0 .- 6.0.*(x_red.*y_red).^2.0 )

    return dy, dx
end