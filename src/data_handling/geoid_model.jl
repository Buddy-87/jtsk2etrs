# https://math.stackexchange.com/questions/4476832/use-fft-in-2-d-space-for-fast-interpolation-of-the-scattered-data


mutable struct Geoid{T<:AbstractFloat}
    name::String
    bmin::T
    bmax::T
    lmin::T
    lmax::T
    db::T
    dl::T
    nrows::Integer
    ncols::Integer
    undolation::Matrix{T}

    function Geoid(name::AbstractString, 
                bmin::T, 
                bmax::T, 
                lmin::T, 
                lmax::T, 
                db::T, 
                dl::T, 
                data::Matrix{T}) where {T<:AbstractFloat}
        this = new{T}()
        
        this.name = name
        this.bmin = bmin
        this.bmax = bmax
        this.lmin = lmin
        this.lmax = lmax
        this.db = db
        this.dl = dl
        
        this.nrows = round(Integer, abs(bmax - bmin)/db)
        this.ncols = round(Integer, abs(lmax - lmin)/dl)

        this.undolation = data

        return this
    end
end # Geoid

function interpolate(g::Geoid{T}, b::T, l::T)::T where (T<:AbstractFloat)
    # check if the coordinates are within the range of the geoid
    if b < g.bmin || b > g.bmax || l < g.lmin || l > g.lmax
        error("Coordinates are outside the range of the geoid.")
    end

    # find the row and column indices of the grid cell containing the coordinates
    row = round(Integer, (b - g.bmin)/g.db)
    col = round(Integer, (l - g.lmin)/g.dl)

    # find the four grid cells surrounding the coordinates
    row_min = max(row - 1, 0)
    row_max = min(row + 1, g.nrows - 1)
    col_min = max(col - 1, 0)
    col_max = min(col + 1, g.ncols - 1)

    # use linear interpolation to find the undulation at the coordinates
    undulation = (g.undolation[row_min, col_min]*(g.bmax - b)*(g.lmax - l) + 
                    g.undolation[row_min, col_max]*(g.bmax - b)*(l - g.lmin) +
                    g.undolation[row_max, col_min]*(b - g.bmin)*(g.lmax - l) +
                    g.undolation[row_max, col_max]*(b - g.bmin)*(l - g.lmin))/
                    (g.bmax - g.bmin)*(g.lmax - g.lmin)

    return undulation
end