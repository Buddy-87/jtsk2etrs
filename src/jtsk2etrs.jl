module jtsk2etrs  
    using ArgParse
    using DelimitedFiles
    using Printf

    # CONSTANTS
    PROGRAM_STATUS_OK = 0
    PROGRAM_STATUS_ERROR = 1
    PROGRAM_ERROR_IN_THE_COMMAND_LINE = 2

    PATH2TABLES = "../data/yx_transformation_tables.hdf5"
    TABLES = ["table_v1005_header", 
              "table_y_3_v1005", 
              "table_x_3_v1005", 
              "table_v1202_header", 
              "table_y_3_v1202", 
              "table_x_3_v1202", 
              "table_v1710_header", 
              "table_y_3_v1710", 
              "table_x_3_v1710"]

    PATH2QGEOIDMODEL = "../data/CR2005.dat"

    VERSION = "0.0.1"

    #
    # cartographic trasnformation
    #
    include("./coordinate_transformation/coordinate_transformation.jl")
    include("./coordinate_transformation/cartographic_projection.jl")
    include("./coordinate_transformation/ellipsoid.jl")
    include("./coordinate_transformation/etrs2jtsk.jl")

    export xyz2blh, blh2xyz, xyz2sph, sph2xyz, xyz2uvw, uvw2xyz
    export bl2jtsk, jtsk2bl, bl2jtsk05, jtsk052bl, interpolate_correction
    export Ellipsoid, grs80, wgs84, bessel
    export get_transformation_parameters, bessel2etrf, etrf2bessel

    #
    # data handling 
    #
    include("./data_handling/geodetic_sphere.jl")
    include("./data_handling/geoid_model.jl")
    include("./data_handling/load_data.jl")
    include("./data_handling/table.jl")

    export geodetic1, geodetic2, length_sphere, lhuilierTheorem, get_right_interval
    export Geoid, interpolate_geoid
    export load_table
    export Table, check_dimensions, bilinear_interpolation, load_raw_table, save_table2hdf5

    #
    # data plot
    #
    include("./plot/plot_raw_data.jl")

    #
    # CLI part
    #
    function convert_etrf2jtsk05( data::Matrix{Any}, table::Table, geoid::Geoid)::Matrix{Any}
        # converting b, l, h  on the wgs84 to cartesian x, y, z coordinates 
        x_etrs, y_etrs, z_etrs = blh2xyz(convert(Vector{Float64}, data[:,2]), 
                                         convert(Vector{Float64}, data[:,3]), 
                                         convert(Vector{Float64}, data[:,4]), wgs84)

        # converting the cartesian x_etrs, y_etrs, z_etrs coordinates of the points related to the wgs84 ellipsoid
        # to the x_bessel, y_bessel, z_bessel related to the bessel ellipsoid
        x_bessel, y_bessel, z_bessel = etrf2bessel(x_etrs, y_etrs, z_etrs)

        # converting the rectangular x_bessel, y_bessel, z_bessel to the geodetic latitude, 
        # longitude and ellipsoidal height using the Bessel ellipsoid, h_bessel is not used 
        # (height conversion is only substracting the geoidal height from ellipsoidal)
        b_bessel, l_bessel, _ = xyz2blh(x_bessel, y_bessel, z_bessel, bessel)

        # converting the geodetic latitude and longitude from (Bessel's ellipsoid) to JTSK-plane
        y_jtsk, x_jtsk = bl2jtsk05(b_bessel, l_bessel)

        n = lastindex(y_jtsk)

        results::Matrix{Any} = Array{Any}(undef, n, 4)

        # converting the ellipsoid height to normal Molodensky's height (h_ell = h_n + zeta)
        # applying the 'tabled' corrections to the final S-JTSK coordinates
        for i=1:n
            dy, dx = bilinear_interpolation(table, y_jtsk[i], x_jtsk[i])
            dzeta = interpolate_geoid(geoid, data[i,2], data[i,3], "sphtriangle")

            results[i,1] = data[i,1]
            results[i,2] = y_jtsk[i] - dy
            results[i,3] = x_jtsk[i] - dx
            results[i,4] = data[i,4] - dzeta
        end

        return results
    end

    function convertjtsk052etrf( data::Matrix{Any}, table::Table, geoid::Geoid)::Matrix{Any}
        n = size(data, 1) 
        results::Matrix{Any} = Array{Any}(undef, n, 4)

        y_jtsk::Vector{Float64} = convert(Vector{Float64}, data[:,2])
        x_jtsk::Vector{Float64} = convert(Vector{Float64}, data[:,3])

        # converting the  normal Molodensky's height (h_n = h_ell - zeta) to ellipsoid height
        # applying the negative 'tabled' corrections to the final S-JTSK coordinates
        for i=1:n
            dy, dx = bilinear_interpolation(table, y_jtsk[i], x_jtsk[i])

            results[i,1] = data[i,1]
            y_jtsk[i] += dy
            x_jtsk[i] += dx            
        end

        # converting y, x JTSK-plane to the geodetic latitude and longitude from (Bessel's ellipsoid) 
        b_bessel, l_bessel = jtsk052bl(y_jtsk, x_jtsk)

        # converting to the cartesian x_etrs, y_etrs, z_etrs coordinates of the points related to the bessel ellipsoid
        x_bessel, y_bessel, z_bessel = blh2xyz(b_bessel, l_bessel, convert(Vector{Float64}, data[:,4]), bessel)

        # using the 7 parameter transformation for converting the x, y, z from
        # Bessel's ellipsoid to the wgs84 (the ETRF reference frame)
        x_etrs, y_etrs, z_etrs = bessel2etrf(x_bessel, y_bessel, z_bessel)

        # conversion to cartesian -> geodetic coordinates
        b_etrs, l_etrs, _ = xyz2blh(x_etrs, y_etrs, z_etrs, wgs84)
        results[:, 2] = b_etrs
        results[:, 3] = l_etrs

        for i=1:n
            dzeta = interpolate_geoid(geoid, b_etrs[i], l_etrs[i], "sphtriangle")
            results[i,4] = data[i,4] + dzeta
        end

        return results
    end

    function save2file(data::Matrix{Any}, fname::AbstractString, computation::AbstractString)
        n = size(data, 1)
        
        open(fname, "w") do f
            for i=1:n
                line = ""
                if computation === "etrf2jtsk"
                    line *= @sprintf "%s    %.3f    %.3f    %.3f\n" data[i,1] data[i,2] data[i,3] data[i,4]
                elseif computation === "jtsk2etrf"
                    line *= @sprintf "%s    %.9F    %.9F    %.3F\n" data[i,1] data[i,2] data[i,3] data[i,4]
                else
                    line *= @sprintf "%s    %.12e    %.12e    %.8e\n" data[i,1] data[i,2] data[i,3] data[i,4]
                end

                write(f, line)
            end
        end
    end


    #
    # command line parser
    # https://argparsejl.readthedocs.io/en/latest/argparse.html
    function parse_commandline()
        settings = ArgParseSettings()
        @add_arg_table! settings begin
            "--file", "-f"
                help = "an option with an argument"
                required = true
            "--transformation", "-t"
                help = "Set the trasnformation either from ETRF -> JTSK or otherwise"
                required = true
            "--output", "-o"
                help = "another option with an argument"
                default = "result.txt"
            "--table"
                help = "Parameter determining which table of 'dy, dx' corrections should be used."
                default = "v1202"
                required = false
            # "arg1"
            #     help = "a positional argument"
            #     required = trues
        end

        #@add_arg_table!(settings, "--help", action => :help, required = false)
        
        return parse_args(settings)
    end


    # ========================== #
    #                            # 
    # C/C++ style MAIN function  #
    #                            #
    # ========================== #
    function main(args::Vector{String})::Integer
        parsed_args = parse_commandline()

        ### arguments and variables for the parse_commandline
        input_coordindates::String = parsed_args["file"] # File with the input coordinates
        transformation::String     = parsed_args["transformation"]
        table::String              = parsed_args["table"]
        output::String             = parsed_args["output"]

        data = load_table(PATH2TABLES, TABLES)
        empty = Array{Float64}(undef, 0,0)
        
        yx_corrections::Table = Table("default", data[1], empty, empty)
        if table === "v1005"
            yx_corrections = Table("table_v1005", data[1], data[2], data[3])
        elseif table === "v1202"
            yx_corrections = Table("table_v1202", data[4], data[5], data[6])
        elseif table === "v1710"
            yx_corrections = Table("table_v1710", data[7], data[8], data[9])
        else
            @warn "Unknown parameter for --table, program ends here."
            return PROGRAM_ERROR_IN_THE_COMMAND_LINE
        end


        if output === "result.txt"
            @warn "The \"--output, -o\" parameter was not provided. The default value is used."
        end

        ### Executable ###
        coordinates::Matrix{Any} = load_points_from_file(input_coordindates)
        qgmodel::Geoid = load_Isgem_model(PATH2QGEOIDMODEL, Float64(1.0))

        results::Matrix{Any} = Matrix{Any}(undef,0,0)
        if transformation === "etrf2jtsk"
            # transformation from the global reference frame into planar coordinate system

            results = convert_etrf2jtsk05(coordinates, yx_corrections, qgmodel)

        elseif transformation === "jtsk2etrf"
            # transformation from the plannar coordinate system into global reference frame

            results = convertjtsk052etrf(coordinates, yx_corrections, qgmodel)
        else
            @warn "Unknown parameter for --transformation, program ends here."
            return PROGRAM_ERROR_IN_THE_COMMAND_LINE
        end

        save2file(results, output, transformation)
    
        return PROGRAM_STATUS_OK
    end


    if abspath(PROGRAM_FILE) == @__FILE__

        # do something only when this file is executed. 
        main(ARGS)

    end
end # module