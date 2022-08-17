"""
Structure to hold the reference ellipsoid parameters.

# Arguments
- `name::AbstractString` the name of the reference ellipsoid
- `a_axis::AbstractFloat` value of the semi-major axis
- `b_axis::AbstractFloat` value of the semi-minor axis

# Examples
```jldoctest
julia> wgs84 = Ellipsoid("World geodetic system 1984", 6378137.000, 6356752.31424518)
```
"""
struct Ellipsoid
    name::String
    a_axis::Float64
    b_axis::Float64
    eccentricity::Float64
end

grs80 = Ellipsoid("Geodetic reference system 1980",6378137.000, 6356752.31424518, 0.006694379990141124)
wgs84 = Ellipsoid("World geodetic system 1984", 6378137.000, 6356752.31424518, 0.006694379990141124)
bessel= Ellipsoid("Bessel", 6377397.155, 6356078.962822, 0.0066743722306109 )