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
                data::Matrix{T})
        this = new()
        
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