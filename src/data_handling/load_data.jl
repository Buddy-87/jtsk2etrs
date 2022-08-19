import .Base.Filesystem.isfile
import DelimitedFiles: readdlm

using HDF5

function file_exist(fname::AbstractString)::Bool
    try        
        open(fname, "r") do _ end
        return true
    catch
        return false
    end
end 

function load_points_from_file(fname::AbstractString)
    if file_exist(fname)
        return readdlm(fname)
    else
        error("Can not access the file '$fname'.") 
    end
end



"""
Load data from a file in HDF5 format. The data stored in this file are:
- `table_v1005_header` : a table v1005 header listing the basic grid info
- `table_y_3_v1005` : a table with the y corrections for the v1005 data
- `table_x_3_v1005` : a table with the x corrections for the v1005 data
- `table_v1202_header` : a table v1202 header listing the basic grid info
- `table_y_3_v1202` : a table with the y corrections for the v1202 data
- `table_x_3_v1202` : a table with the x corrections for the v1202 data
- `table_v1710_header` : a table v1710 header listing the basic grid info
- `table_y_3_v1710` : a table with the y corrections for the v1710 data
- `table_x_3_v1710` : a table with the x corrections for the v1710 data
- #TODO geoid header vector
- `CZE_qgeoid_2005_v1005`` : a table with the undulation corrections for the CZE quasigeoid 2005 model

# Arguments
- `fname`: name of the file to load
# Returns
- `data`: a dictionary with the data

# Example
```
data = load_table("yx_transformation_tables.hdf5", ["table_v1005_header", "table_y_3_v1005", "table_x_3_v1005"])
```
"""
function load_table(fname::AbstractString, data_req::AbstractVector)
    if isempty(data_req)
        # if data_req is empty returns empty Vector{Any}
        return Any[]
    end

    if file_exist(fname)
        data = Any[]

        h5open(fname, "r") do file
            for read_data in data_req
                if read_data == "table_v1005_header"
                    push!(data, read(file, "table_v1005_header"))
                elseif read_data == "table_y_3_v1005"
                    push!(data, read(file, "table_y_3_v1005"))
                elseif read_data == "table_x_3_v1005"
                    push!(data, read(file, "table_x_3_v1005"))
                elseif read_data == "table_v1202_header"
                    push!(data, read(file, "table_v1202_header"))
                elseif read_data == "table_y_3_v1202"
                    push!(data, read(file, "table_y_3_v1202"))
                elseif read_data == "table_x_3_v1202"
                    push!(data, read(file, "table_x_3_v1202"))
                elseif read_data == "table_v1710_header"
                    push!(data, read(file, "table_v1710_header"))
                elseif read_data == "table_y_3_v1710"
                    push!(data, read(file, "table_y_3_v1710"))
                elseif read_data == "table_x_3_v1710"
                    push!(data, read(file, "table_x_3_v1710"))
                elseif read_data == "CZE_qgeoid_2005_v1005"
                    push!(data, read(file, "CZE_qgeoid_2005_v1005"))
                else
                    @warn "Requested \"$read_data\" are not supported. See the documentation."
                end
            end
        end

        return data
    else 
        error("Can not access the file '$fname'.") 
    end
end 