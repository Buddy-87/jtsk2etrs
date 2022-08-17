using DelimitedFiles
using HDF5

include("./load_data.jl")

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

        nr_y, nc_y = size(yCorrections)
        nr_x, nc_x = size(xCorrections)

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
Load raw table of x-corrections and y-corrections as they are provided by the Czech 
Cadastral office.

# Arguments

# Returns

# Example
"""
function load_raw_table(fname::AbstractString, modelname::AbstractString = "unknown")::Table
    if file_exist(fname)
        table = readdlm(fname)

        nr = size(table)[1]

        ymin::Float64 = minimum(table[:,1]) 
        ymax::Float64 = maximum(table[:,1])
        xmin::Float64 = minimum(table[:,2])
        xmax::Float64 = maximum(table[:,2])

        dy::Float64 = 2000.0 # constant value for tables table_yx_3_v**** provided by the CUZK
        dx::Float64 = 2000.0 

        table_nr::Integer = round(Integer, abs(xmax - xmin)/dx) + 1 
        table_nc::Integer = round(Integer, abs(ymax - ymin)/dy) + 1 

        y_table = Array{Float64}(undef, table_nr, table_nc)
        x_table = Array{Float64}(undef, table_nr, table_nc)

        for r = 1:nr
            row::Integer = round(Integer, (table[r,2] - xmin)/dx) + 1
            col::Integer = round(Integer, (table[r,1] - ymin)/dy) + 1

            y_table[row, col] = table[r,3]
            x_table[row, col] = table[r,4]
        end

        header = [ymin, ymax, xmin, xmax, dy, dx]

        t::Table{Float64} = Table(modelname, header, y_table, x_table)
        return t
    else
        error("File $fname does not exist.")
    end
end

"""
Saves/appends the table to hdf file 

# Arguments

# Returns
- `true` if the process is successful, `false` otherwise
# Example
```
# save to file 'yx_transformation_tables.HDF5' data stored in 'table_v1005'
julia> save_table2hdf5("yx_transformation_tables.HDF5", table_v1005, ["table_v1005_header"; "table_y_3_v1005"; "table_x_3_v1005"])
true
```

"""
function save_table2hdf5(outname::AbstractString, table::Table, dataNames::AbstractVector)::Bool
    header = [table.ymin, table.ymax, table.xmin, table.xmax, table.dy, table.dx]

    try
        h5open(outname ,isfile(outname) ? "r+" : "w") do file
            write(file, dataNames[1], header)
            write(file, dataNames[2], table.yCorrections)
            write(file, dataNames[3], table.xCorrections)
        end

        return true
    catch
        return false
    end
end


"""
# Arguments
"""
function check_dimensions(this::Table{T}, data::Matrix{T}) where {T<:AbstractFloat}
    nr_y, nc_y = this.yCorrections.shape
    nr_x, nc_x = this.xCorrections.shape

    if nr_y != nr_x || nc_y != nc_x
        error("yCorrections and xCorrections must have the same shape")
    end

    if nr_y != data.shape[1] || nc_y != data.shape[2]
        error("data must have the same shape as yCorrections and xCorrections")
    end
end

"""
In mathematics, bilinear interpolation is a method for interpolating functions of two variables
(e.g., x and y) using repeated linear interpolation. It is usually applied to functions sampled 
on a 2D rectilinear grid, though it can be generalized to functions defined on the vertices
of (a mesh of) arbitrary convex quadrilaterals.

https://en.wikipedia.org/wiki/Bilinear_interpolation

# Arguments
Suppose that we want to find the value of the unknown function f at the point (x, y). It is assumed that we know the value of f at the four points f11 = (x1, y1), f12 = (x1, y2), f21 = (x2, y1), and f22 = (x2, y2).

- `x::T`: x-coordinate 
- `y::T`: y-coordinate 
- `x1::T`: x-coordinate 
- `x2::T`: x-coordinate 
- `y1::T`: y-coordinate 
- `y2::T`: y-coordinate 
- `f11::T`: function value in the node f11 = f(x1, y1) 
- `f12::T`: function value in the node f12 = f(x1, y2)  
- `f21::T`: function value in the node f21 = f(x2, y1)  
- `f22::T`: function value in the node f22 = f(x2, y2)  

# Returns
- Interpolated value

# Example
```
julia> x1, x2, y1, y2 = 14.0, 15.0, 21.0, 20.0
julia> f11, f12, f21, f22 = 162.0, 91.0, 95.0, 210.0
bilinear_interpolation(14.5, 20.2, x1, x2, y1, y2, f11, f12, f21, f22)
146.09999999999854
```
"""
function bilinear_interpolation(x::T, y::T, x1::T, x2::T, y1::T, y2::T, f11::T, f12::T, f21::T, f22::T)::T where {T<:AbstractFloat}
    d0::T = (x1-x2)*(y1-y2);

    a00::T =  f11*x2*y2 - f12*x2*y1 - f21*x1*y2 + f22*x1*y1;
    a10::T = -f11*y2    + f12*y1    + f21*y2    - f22*y1;   
    a01::T = -f11*x2    + f12*x2    + f21*x1    - f22*x1;
    a11::T =  f11       - f12       - f21       + f22;

    return (a00 + a10*x + a01*y + a11*x*y)/d0;
end

"""
Interpolates the value from the 'Table' containing the corrections that are related to 
'y' and 'x' coordinates. The function checks if the requested coordinates are in the
bounds of the 'Table'.
# Arguments
Suppose that we want to find the value of the unknown function f at the point (x, y). It is assumed that we know the value of f at the four points f11 = (x1, y1), f12 = (x1, y2), f21 = (x2, y1), and f22 = (x2, y2).

- `table::Table` table with the corrections @see Table
- `x::T`: x-coordinate 
- `y::T`: y-coordinate 

# Returns
- Interpolated values in pair `(iy, ix)`, returns `(0.0, 0.0)` if the requested Coordinates
are outside the table bounds.

# Example
```

julia> data1 = load_table("yx_transformation_tables.hdf5", ["table_v1005_header", "table_y_3_v1005", "table_x_3_v1005"])
julia> t1 = Table("table_v1005", data1[1], data1[2], data1[3])
julia> bilinear_interpolation(t1, 512464.123, 1089123.654)
(0.15616191491377354, 0.11497914209079742)

```
"""
function bilinear_interpolation(table::Table, y::T, x::T) where {T<:AbstractFloat}
    if y >=table.ymin && y <= table.ymax && x >= table.xmin && x <= table.xmax
        row::Integer = round(Integer, (x - table.xmin)/table.dx) + 1
        col::Integer = round(Integer, (y - table.ymin)/table.dy) + 1

        x1::T = table.xmin + T(row-1) * table.dx
        x2::T = x1 + table.dx

        y1::T = table.ymin + T(col-1) * table.dy
        y2::T = y1 + table.dy

        fy11::T = table.yCorrections[row, col]
        fy12::T = table.yCorrections[row, col+1]
        fy21::T = table.yCorrections[row+1, col]
        fy22::T = table.yCorrections[row+1, col+1]

        fx11::T = table.xCorrections[row, col]
        fx12::T = table.xCorrections[row, col+1]
        fx21::T = table.xCorrections[row+1, col]
        fx22::T = table.xCorrections[row+1, col+1]

        y_int::T = bilinear_interpolation(x, y, x1, x2, y1, y2, fy11, fy12, fy21, fy22)
        x_int::T = bilinear_interpolation(x, y, x1, x2, y1, y2, fx11, fx12, fx21, fx22)

        return y_int, x_int
    else
        return 0.0, 0.0
    end
end