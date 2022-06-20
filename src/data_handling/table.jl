mutable struct Table{T<:AbstractFloat}
    name::String
    ymin::T
    ymax::T
    xmin::T
    xmax::T
    dy::T
    dx::T
    nrows::Integer
    ncols::Integer
    yCorrections::Matrix{T}
    xCorrections::Matrix{T}

    function Table(name::AbstractString, 
                   header::Vector{T},  
                   yCorrections::Matrix{T},
                   xCorrections::Matrix{T}) where {T<:AbstractFloat}
        this = new{T}()
        
        this.name = name

        this.ymin = header[1]
        this.ymax = header[2]
        this.xmin = header[3]
        this.xmax = header[4]
        this.dy = header[5]
        this.dx = header[6]

        nr_y, nc_y = yCorrections.shape
        nr_x, nc_x = xCorrections.shape

        if nr_y != nr_x || nc_y != nc_x
            error("yCorrections and xCorrections must have the same shape")
        end

        this.nrows = nr_y
        this.ncols = nc_y

        this.yCorrections = yCorrections
        this.xCorrections = xCorrections        

        return this
    end
end

"""
# Arguments
"""
function check_dimensions(this::Table{T}, data::Matrix{T}) where (T<:AbstractFloat)
    nr_y, nc_y = this.yCorrections.shape
    nr_x, nc_x = this.xCorrections.shape

    if nr_y != nr_x || nc_y != nc_x
        error("yCorrections and xCorrections must have the same shape")
    end

    if nr_y != data.shape[1] || nc_y != data.shape[2]
        error("data must have the same shape as yCorrections and xCorrections")
    end
end