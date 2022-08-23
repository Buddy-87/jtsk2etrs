# https://math.stackexchange.com/questions/4476832/use-fft-in-2-d-space-for-fast-interpolation-of-the-scattered-data

using DelimitedFiles

import LinearAlgebra: det

include("./load_data.jl")
include("./geodetic_sphere.jl")

function load_Isgem_model(fname::AbstractString, numerical_type::T)::Geoid{T} where {T<:AbstractFloat}
    name::String      = "" # name of the model
    type::String      = "" # type of the model geoid/quasi-geoid
    units::String     = "" # units like [m]
    reference::String = "" # reference ellipsoid

    lat_min::T = -90.0  # minimum value of the latitude
    lat_max::T = +90.0  # maximum value of the latitude
    lon_min::T = +0.0   # minimum value of the longitude
    lon_max::T = +360.0 # maximum value of the longitude
    dlat::T = 1.0       # difference between the latitudal lines, e.g. step
    dlon::T = 1.0       # difference between the longitudal lines, e.g. step
    nrows::Integer = 0  # number of nrows
    ncols::Integer = 0  # number of columns
    nodata::T = -9999   # difference between the longitudal lines, e.g. step

    eoh::Integer = 1 # number of lines to reach end of header
    if file_exist(fname)
        open(fname) do io
            while !eof(io)
                line::String = readline(io)
                # println("\"$line\"")

                words = split(line)
                
                if occursin("model name", line) 
                    name = last(words)
                elseif occursin("model type", line) 
                    type = last(words)
                elseif occursin("units", line)      
                    units = last(words)
                elseif occursin("reference", line)  
                    reference = last(words)
                elseif occursin("lat min", line)    
                    lat_min = parse(T, last(words))
                elseif occursin("lat max", line)    
                    lat_max= parse(T, last(words))
                elseif occursin("lon min", line)    
                    lon_min = parse(T, last(words))
                elseif occursin("lon max", line)    
                    lon_max = parse(T, last(words))
                elseif occursin("delta lat", line)  
                    dlat = parse(T, last(words))
                elseif occursin("delta lon", line)  
                    dlon = parse(T, last(words))
                elseif occursin("nrows", line)      
                    nrows = parse(Int64, last(words))
                elseif occursin("ncols", line)      
                    ncols = parse(Int64, last(words))
                elseif occursin("nodata", line)     
                    nodata = parse(T, last(words))
                end
 
                # end of the header
                occursin("end_of_head", line) && break
                eoh += 1
            end
        end

        data::Matrix{T} = readdlm(fname, skipstart=eoh)
        nr, nc = size(data)
        
        if nr != nrows && nc != ncols
            error("The header info does not correspond with the actual data. Invalid matrix size!")
        end

        geoid_model = Geoid(name, 
                            lat_min,
                            lat_max,
                            lon_min,
                            lon_max,
                            dlat,
                            dlon,
                            nodata,
                            flipud(data))

        return geoid_model
    else
        error("File $fname does not exist.")
    end
end

mutable struct Geoid{T<:AbstractFloat}
    name::String
    bmin::T
    bmax::T
    lmin::T
    lmax::T
    db::T
    dl::T
    nodata::T
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
                nodata::T,
                data::Matrix{T}) where {T<:AbstractFloat}
        this = new{T}()
        
        this.name = name
        this.bmin = bmin
        this.bmax = bmax
        this.lmin = lmin
        this.lmax = lmax
        this.db = db
        this.dl = dl
        this.nodata = nodata
        
        this.nrows = round(Integer, abs(bmax - bmin)/db)
        this.ncols = round(Integer, abs(lmax - lmin)/dl)

        this.undolation = data

        return this
    end
end # Geoid

function interpolate_geoid(g::Geoid{T}, b::T, l::T, algorithm::AbstractString="nearest") where (T<:AbstractFloat)
    # check if the coordinates are within the range of the geoid
    if b < g.bmin || b > g.bmax || l < g.lmin || l > g.lmax
        error("Coordinates are outside the range of the geoid.")
    end

    # find the row and column indices of the grid cell containing the coordinates
    row = round(Integer, (b - g.bmin)/g.db)
    col = round(Integer, (l - g.lmin)/g.dl)

    # find the four grid cells surrounding the coordinates
    row_min = max(row - 1, 0)
    row_max = min(row, g.nrows)
    col_min = max(col - 1, 0)
    col_max = min(col, g.ncols)

    undulation::T = 0.0;
    bcen::T = T(row * g.db + g.bmin)
    lcen::T = T(col * g.dl + g.lmin)
    bdif::T = b - (bcen + g.db/2.0);
    ldif::T = l - (lcen + g.db/2.0);

    phi_list = Vector{T}(undef, 4)
    lam_list = Vector{T}(undef, 4)
    fvalues = Vector{T}(undef, 4)


    result::T = 0.0; w::T = 0.0; pw::T = 0.0;
    if algorithm === "sphtriangle"
        #= Linear interpolation =#
        if ( bdif < 0 && ldif < 0 ) 
            # bottom left, bottom right, upper left (bl, br, ul)
            phi_list[1] = bcen - g.db/2.0;
            phi_list[2] = bcen - g.db/2.0;
            phi_list[3] = bcen + g.db/2.0;
            phi_list[4] = bcen - g.db/2.0;

            lam_list[1] = lcen - g.dl/2.0;
            lam_list[2] = lcen + g.dl/2.0;
            lam_list[3] = lcen - g.dl/2.0;
            lam_list[4] = lcen - g.dl/2.0;

            fvalues[1] = g.undolation[row+1, col]
            fvalues[2] = g.undolation[row, col]
            fvalues[3] = g.undolation[row, col+1]

        elseif (bdif >= 0 && ldif < 0 )
            # bl, ul, ur
            phi_list[1] = bcen - g.db/2.0;
            phi_list[2] = bcen + g.db/2.0;
            phi_list[3] = bcen + g.db/2.0;
            phi_list[4] = bcen - g.db/2.0;

            lam_list[1] = lcen - g.dl/2.0;
            lam_list[2] = lcen - g.dl/2.0;
            lam_list[3] = lcen + g.dl/2.0;
            lam_list[4] = lcen - g.dl/2.0;

            fvalues[1] = g.undolation[row+1, col+1]
            fvalues[2] = g.undolation[row, col]
            fvalues[3] = g.undolation[row+1, col]
        elseif ( bdif >= 0 && ldif >= 0 )
            # br, ul, ur
            phi_list[1] = bcen - g.db/2.0;
            phi_list[2] = bcen + g.db/2.0;
            phi_list[3] = bcen + g.db/2.0;
            phi_list[4] = bcen - g.db/2.0;

            lam_list[1] = lcen + g.dl/2.0;
            lam_list[2] = lcen - g.dl/2.0;
            lam_list[3] = lcen + g.dl/2.0;
            lam_list[4] = lcen + g.dl/2.0;

            fvalues[1] = g.undolation[row+1, col+1];
            fvalues[2] = g.undolation[row, col+1];
            fvalues[3] = g.undolation[row+1, col];
        else
            # bl, br, ur
            phi_list[1] = bcen - g.db/2.0;
            phi_list[2] = bcen - g.db/2.0;
            phi_list[3] = bcen + g.db/2.0;
            phi_list[4] = bcen - g.db/2.0;

            lam_list[1] = lcen - g.dl/2.0;
            lam_list[2] = lcen + g.dl/2.0;
            lam_list[3] = lcen + g.dl/2.0;
            lam_list[4] = lcen - g.dl/2.0;

            fvalues[1] =  g.undolation[row+1, col+1]
            fvalues[2] =  g.undolation[row, col]
            fvalues[3] =  g.undolation[row, col+1]
        end

        s1::T = length_sphere( phi_list[1], lam_list[1],  phi_list[2], lam_list[2]);
        s2::T = length_sphere( phi_list[2], lam_list[2],  phi_list[3], lam_list[3]);
        s3::T = length_sphere( phi_list[3], lam_list[3],  phi_list[1], lam_list[1]);

        surface_triangle::T = lhuilierTheorem(s1, s2, s3);

        for i=1:3
            s1 = length_sphere( phi_list[i], lam_list[i],  phi_list[i+1], lam_list[i+1]);
            s2 = length_sphere( phi_list[i], lam_list[i],  b, l);
            s3 = length_sphere( b, l,  phi_list[i+1], lam_list[i+1]);

            surface::T = lhuilierTheorem(s1, s2, s3); # in radians, of given triangle
            weight::T = surface/surface_triangle;

            w  += weight;
            pw += weight * fvalues[i];
        end

        result = pw/w;

        if (isnan(result) || isinf(result))
            return g.nodata
        end

        return result
    else # algorithm === "nearest"
        #=  Distance-based interpolation =#
        phi_list[1] = bcen - g.db/2.0;
        phi_list[2] = bcen - g.db/2.0;
        phi_list[3] = bcen + g.db/2.0;
        phi_list[4] = bcen + g.db/2.0;

        lam_list[1] = lcen - g.dl/2.0;
        lam_list[2] = lcen + g.dl/2.0;
        lam_list[3] = lcen - g.dl/2.0;
        lam_list[4] = lcen + g.dl/2.0;

        fvalues[1] = g.undolation[row, col]
        fvalues[2] = g.undolation[row, col+1]
        fvalues[3] = g.undolation[row+1, col]
        fvalues[4] = g.undolation[row+1, col+1]

        for i=1:4
            si::T = length_sphere(b,l, phi_list[i], lam_list[i])
            wi::T = 1.0/si^2.0

            w  += wi;
            pw += wi * fvalues[i];
        end
        
        result = pw/w;

        if (isnan(result) || isinf(result))
            return g.nodata
        end

        return result
    end

    return undulation
end

"""
Creates a deep copy of the matrix and flips it upside down
"""
function flipud(data::Matrix{T}) where {T<:AbstractFloat}
    nr, nc = size(data)

    dataud = Array{T}(undef, nr, nc)

    for i=1:nr
        dataud[i, :] = deepcopy(data[end-(i-1), :])
    end

    return dataud
end