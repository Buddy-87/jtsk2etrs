module ellipsoid
    export Ellipsoid
    export wgs84
    export grs80

    """
    Structure to hold the reference ellipsoid parameters.

    # Arguments
    - `name::AbstractString` the name of the reference ellipsoid
    - `a_axis::Real` value of the semi-major axis
    - `b_axis::Real` value of the semi-minor axis

    # Examples
    ```jldoctest
    julia> Ellipsoid("World geodetic system 1984", 6378137.000, 6356752.31424518)
    ```
    """
    mutable struct Ellipsoid
        name::String
        a_axis::Float64
        b_axis::Float64
        eccentricity::Float64

        function Ellipsoid(name::AbstractString, a_axis::Real, b_axis::Real)
            this = new()

            this.name = name
            this.a_axis = a_axis
            this.b_axis = b_axis

            this.eccentricity = (a_axis^2 - b_axis^2)/a_axis^2

            return this
        end
    end
   
    
    const wgs84 = Ellipsoid("World geodetic system 1984", 6378137.000, 6356752.31424518)
    const grs80 = Ellipsoid("Geodetic reference system 1980",6378137.000, 6356752.31424518)
end # ellipsoid module