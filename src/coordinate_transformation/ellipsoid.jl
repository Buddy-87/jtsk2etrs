module ellipsoid
    export Ellipsoid
    export wgs84
    export grs80

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
mutable struct Ellipsoid{T<:AbstractFloat}
    name::String
    a_axis::T
    b_axis::T
    eccentricity::T

    function Ellipsoid(name::AbstractString, a_axis::T, b_axis::T) where {T<:AbstractFloat}
        this = new{T}()

        this.name = name
        this.a_axis = a_axis
        this.b_axis = b_axis

        this.eccentricity = (a_axis^2 - b_axis^2)/a_axis^2

        return this
    end
end


grs80 = Ellipsoid("Geodetic reference system 1980",6378137.000, 6356752.31424518)
wgs84 = Ellipsoid("World geodetic system 1984", 6378137.000, 6356752.31424518)
bessel= Ellipsoid("Bessel", 6377397.155,  6356078.962822 )

end # ellipsoid module