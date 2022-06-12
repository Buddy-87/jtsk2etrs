"""
Return the transformation parameter for transformation between the
ETRF2000 reference ellipsoid and the Bessel reference ellipsoid. 
The parameter are usable for the function `etrf2bessel` and function
`bessel2etrf`.

The parameters are different for the inverse transformation (from Bessel to ETRF2000).
This is mainly caused due to the relativelly high angular distance between the
the  ETRF2000 and S-JTKS (Bessel).

# Arguments
- `transformation::Integer` transformation type, 1 for forward transformation,
  -1 for inverse transformation
# Returns
- `mscale::T` scale factor of the transformation
- `rmat::Matrix{T}` rotation matrix
- `tmat::Vector{T}` translation vector
"""
function get_transformation_parameters::T(transformation::Integer) where (T<:Real)
    # rho is constant for converting the arcseconds to radians
    rho::T = 206264.80624709636;

    if transformation == 1
        # Parameters for converting ETRF2000 to JTSK05
        t1::T = -572.203
        t2::T = -85.328 
        t3::T = -461.934 

        p1::T =  5.24832714/rho
        p2::T =  1.52900087/rho
        p3::T =  4.97311727/rho

        mscale::T = 1+(-3.5393e-6)
        rmat::Matrix{T} = [1 p1 -p2; -p1 1 p3; p2 -p3 1];
        tmat::Matrix{T} = [t1; t2; t3];

        return mscale, rmat, tmat
    elseif transformation == -1
        # Parameters for converting JTSK05 to ETRF2000
        t1::T = 572.213
        t2::T = 85.334
        t3::T = 461.940
        p1::T = -5.24836073/rho
        p2::T = -1.52899176/rho
        p3::T = -4.97316164/rho

        mscale::T = 1+(3.5378e-6)
        rmat::Matrix{T} = [1 p1 -p2; -p1 1 p3; p2 -p3 1];
        tmat::Matrix{T} = [t1; t2; t3];

        return mscale, rmat, tmat
    else
        error("Transformation type not supported")
    end
end

"""
Converting the cartesian coordinates x, y, z form the Bessel
ellipsoid to the ETRF2000. Using the 7-parameter transformation.
The parameters are: 
- m for the scale factor
- p1, p2, p3 are related parameters to the rotation matrix
- t1, t2, t3 are three parameters of the translation vector

# Arguments
- `x::T` x-coordinate in Bessel reference frame
- `y::T` y-coordinate in Bessel reference frame
- `z::T` z-coordinate in Bessel reference frame
# Returns
- `x_etrf::T` x-coordinate in ETRF2000 reference frame
- `y_etrf::T` y-coordinate in ETRF2000 reference frame
- `z_etrf::T` z-coordinate in ETRF2000 reference frame
# Examples
```jldoctest
julia> 
```
"""
function bessel2etrf::T( x::T, y::T, z::T ) where (T<:Real)
    mscale, rmat, tmat = get_transformation_parameters(1);

    # Transformation
    x_etrf::T = mscale*(rmat[1,1]*x + rmat[1,2]*y + rmat[1,3]*z + tmat[1])
    y_etrf::T = mscale*(rmat[2,1]*x + rmat[2,2]*y + rmat[2,3]*z + tmat[2])
    z_etrf::T = mscale*(rmat[3,1]*x + rmat[3,2]*y + rmat[3,3]*z + tmat[3])

    return x_etrf, y_etrf, z_etrf
end

"""
Converting the cartesian coordinates x, y, z form the ETRF2000
ellipsoid to the Bessel ellipsoid. Using the 7-parameter transformation.
The parameters are:
- m for the scale factor
- p1, p2, p3 are related parameters to the rotation matrix
- t1, t2, t3 are three parameters of the translation vector
# Arguments
- `x::T` x-coordinate in ETRF2000 reference frame
- `y::T` y-coordinate in ETRF2000 reference frame
- `z::T` z-coordinate in ETRF2000 reference frame
# Returns
- `x_bessel::T` x-coordinate in Bessel reference frame
- `y_bessel::T` y-coordinate in Bessel reference frame
- `z_bessel::T` z-coordinate in Bessel reference frame
# Examples
```jldoctest
julia>
```
"""
function etrf2bessel::T() where (T<:Real)
    mscale, rmat, tmat = get_transformation_parameters(-1);

    # Transformation
    x_bessel::T = mscale*(rmat[1,1]*x + rmat[1,2]*y + rmat[1,3]*z + tmat[1])
    y_bessel::T = mscale*(rmat[2,1]*x + rmat[2,2]*y + rmat[2,3]*z + tmat[2])
    z_bessel::T = mscale*(rmat[3,1]*x + rmat[3,2]*y + rmat[3,3]*z + tmat[3])

    return x_bessel, y_bessel, z_bessel
    
end

"""
Converting the cartesian coordinates x, y, z form the Bessel ellipsoid to
the ETRF2000 ellipsoid. Using the 7-parameter transformation.
The parameters are:
- m for the scale factor
- p1, p2, p3 are related parameters to the rotation matrix
- t1, t2, t3 are three parameters of the translation vector
# Arguments
- `x::Vector{T}` x-coordinates in column vector in Bessel reference frame
- `y::Vector{T}` y-coordinates in column vector in Bessel reference frame
- `z::Vector{T}` z-coordinates in column vector in Bessel reference frame
# Returns
- `x_etrf::Vector{T}` x-coordinates in column vector in ETRF2000 reference frame
- `y_etrf::Vector{T}` y-coordinates in column vector in ETRF2000 reference frame
- `z_etrf::Vector{T}` z-coordinates in column vector in ETRF2000 reference frame
# Examples
```jldoctest
julia>
```
"""
function bessel2etrf::T( x::Vector{T}, y::Vector{T}, z::Vector{T} ) where (T<:Real)
    # create a matrix with the cartesian coordinates
    xyz::Matrix{T} = [x, y, z]

    # transformation
    mscale, rmat, tmat = get_transformation_parameters(1);

    xyz_etrf = transpose(mscale*(rmat* transpose(xyz)))
    (x_etrf, y_etrf, z_etrf) = [xyz_etrf[:,x] for x in 1:size(xyz_etrf,1)]

    return x_etrf + tmat[1], y_etrf + tmat[2], z_etrf + tmat[3]
end

"""
Converting the cartesian coordinates x, y, z form the ETRF2000 ellipsoid to
the Bessel ellipsoid. Using the 7-parameter transformation.
The parameters are:
- m for the scale factor
- p1, p2, p3 are related parameters to the rotation matrix
- t1, t2, t3 are three parameters of the translation vector
# Arguments
- `x::Vector{T}` x-coordinates in column vector in ETRF2000 reference frame
- `y::Vector{T}` y-coordinates in column vector in ETRF2000 reference frame
- `z::Vector{T}` z-coordinates in column vector in ETRF2000 reference frame
# Returns
- `x_bessel::Vector{T}` x-coordinates in column vector in Bessel reference frame
- `y_bessel::Vector{T}` y-coordinates in column vector in Bessel reference frame
- `z_bessel::Vector{T}` z-coordinates in column vector in Bessel reference frame
# Examples
```jldoctest
julia>
```
"""
function etrf2bessel::T( x::Vector{T}, y::Vector{T}, z::Vector{T} ) where (T<:Real)
    # create a matrix with the cartesian coordinates
    xyz::Matrix{T} = [x, y, z]

    # transformation
    mscale, rmat, tmat = get_transformation_parameters(-1);

    xyz_bessel = transpose(mscale*(rmat* transpose(xyz)))
    (x_bessel, y_bessel, z_bessel) = [xyz_bessel[:,x] for x in 1:size(xyz_bessel,1)]

    return x_bessel + tmat[1], y_bessel + tmat[2], z_bessel + tmat[3]
end