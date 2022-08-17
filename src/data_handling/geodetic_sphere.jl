import .Base: sin, cos, tan, asin, deg2rad, rad2deg, abs
             acos, atan, sqrt, sign, isinf, isnan

function geodetic1(phi1::T, lambda1::T, alpha1::T, s1::T, radius::T = 1.0) where {T<:AbstractFloat}
    b1::T = deg2rad(phi1);
    l1::T = deg2rad(lambda1);
    a1::T = deg2rad(alpha1);
    sr::T = s1/radius;

    b2::T = asin(sin(b1)*cos(sr)+ cos(b1)*sin(sr)*cos(a1));
    deltal::T = sign( sin(a1) )* abs(acos((cos(sr)-sin(b1)*sin(b2))/(cos(b1)*cos(b2))));

    a::T = (cos(b1)*sin(a1)/(cos(b2)));
    b::T = -1.0*cos(a1)*cos(deltal) + sin(a1)*sin(deltal)*sin(b1);

    a2::T = 2.0*pi - atan(a, b);

    a1 = get_right_interval( a1 );
    a2 = get_right_interval( a2 );

    phi2::T = rad2dms( b2 );
    lambda2::T = rad2dms( l1+deltal );
    alpha2::T = rad2dms( a2 );

    return phi2, lambda2, alpha2
end

function geodetic2(phi1::T, lambda1::T, phi2::T, lambda2::T, radius::T = 1.0) where {T<:AbstractFloat}
    b1::T = deg2rad(phi1)
    l1::T = deg2rad(lambda1)
    b2::T = deg2rad(phi2)
    l2::T = deg2rad(lambda2)

    dlambda::T = l2 - l1
    s1::T = acos(sin(b1)*sin(b2)+cos(b1)*cos(b2)*cos(dlambda))
    a::T =  atan(cos((b1-b2)/2.),sin((b2+b1)/2.)*tan(dlambda/2))
    b::T =  atan(sin((b1-b2)/2.),cos((b2+b1)/2.)*tan(dlambda/2))

    a1::T = (a+b)
    a2::T = 2.0*pi- (a-b)
    # Conditions azimuth in <0,2 pi>
    a1 = get_right_interval( a1 )
    a2 = get_right_interval( a2 )

    # results
    return s1*radius, rad2deg(a1), rad2deg(a2)
end

function length_sphere(phi1::T, lam1::T, phi2::T, lam2::T, radius::T = 1.0) where {T<:AbstractFloat}
    b1::T = deg2rad(phi1)
    b2::T = deg2rad(phi2)
    dlambda::T = deg2rad(lam2-lam1)

    return radius*acos(sin(b1)*sin(b2)+cos(b1)*cos(b2)*cos(dlambda));
end

function lhuilierTheorem(a::T, b::T, c::T, rsphere::T = 1.0) where {T<:AbstractFloat}
#=  a,b,c lengths of the triangle sides
    r radius of a reference sphere =#
    a /= rsphere;
    b /= rsphere;
    c /= rsphere;

    s::T = .5*(a+b+c); # semiperimeter
    param::T = tan(s/2.0)*tan((s-a)/2.0)* tan((s-b)/2.0)* tan((s-c)/2.0)

    if param < 0 
        return 0.0
    end

    return 4.0*rsphere*rsphere*atan(sqrt(param)); # in radias
end

function get_right_interval(x::T) where {T<:AbstractFloat}
    if x < 0 
        return x + 2.0*pi
    elseif x > 2.0*pi
        return x - 2.0 * pi
    end

    return x
end