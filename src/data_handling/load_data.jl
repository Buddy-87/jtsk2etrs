"""
Load data from a file in HDF5 format. The data stored in this file are:
- `Data_names` : a list of strings with the names of the data
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
"""
function load_table(fname::AbstractString)
    h5open(fname, "r") do file
        vector_names = read(file, "Data_names")
        table_v1202_header = read(file, "table_v1202_header")
        table_y_3_v1202 = read(file, "table_y_3_v1202")
        table_x_3_v1202 = read(file, "table_x_3_v1202")
        table_v1710_header = read(file, "table_v1710_header")
        table_y_3_v1710 = read(file, "table_y_3_v1710")
        table_x_3_v1710 = read(file, "table_x_3_v1710")
        CZE_qgeoid_2005_v1005 = read(file, "CZE_qgeoid_2005_v1005")
    end
    return vector_names, table_v1202_header, table_y_3_v1202, table_x_3_v1202, table_v1710_header, table_y_3_v1710, table_x_3_v1710, CZE_qgeoid_2005_v1005
end