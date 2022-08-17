module jtsk2etrs    
    using ArgParse

    # constants
    PROGRAM_STATUS_OK = 0
    PROGRAM_STATUS_ERROR = 1

    #
    # cartographic trasnformation
    #
    include("coordinate_transformation/coordinate_transformation.jl")
    include("coordinate_transformation/cartographic_projection.jl")
    include("coordinate_transformation/ellipsoid.jl")
    include("coordinate_transformation/etrs2jtsk.jl")

    export xyz2blh, blh2xyz, xyz2sph, sph2xyz, xyz2uvw, uvw2xyz
    export bl2jtsk, jtsk2bl, bl2jtsk05, jtsk052bl, interpolate_correction
    export Ellipsoid, grs80, wgs84, bessel
    export get_transformation_parameters, bessel2etrf, etrf2bessel

    #
    # data handling 
    #
    include("data_handling/geodetic_sphere.jl")
    include("data_handling/geoid_model.jl")
    include("data_handling/load_data.jl")
    include("data_handling/table.jl")

    export geodetic1, geodetic2, length_sphere, lhuilierTheorem, get_right_interval
    export Geoid, interpolate
    export load_table
    export Table, check_dimensions, bilinear_interpolation, load_raw_table, save_table2hdf5

    #
    # data plot
    #
    include("plot/plot_raw_data.jl")

    #
    # command line parser
    # https://argparsejl.readthedocs.io/en/latest/argparse.html
    function parse_commandline()
        settings = ArgParseSettings()
        @add_arg_table settings begin
            "--file", "-f"
                help = "an option with an argument"
            "--trasnformation", "-t"
                help = "Set the trasnformation either from ETRF -> JTSK or otherwise"
            "--output", "-o"
                help = "another option with an argument"
                arg_type = String
                default = "result.txt"
            # "arg1"
            #     help = "a positional argument"
            #     required = true
        end

        #@add_arg_table!(settings, "--help", action => :help, required = false)
        
        return parse_args(["--help", "--version"], settings)
    end


    function main(args::Vector{String})::Integer
        println("Executing the program \"jtsk2etrs\".")
        parsed_args = parse_commandline()
        for (arg,val) in parsed_args
            print("  $arg  =>  ")
            show(val)
            println()
        end
        ### create a command line parser ###

        


        return PROGRAM_STATUS_OK
    end


    if abspath(PROGRAM_FILE) == @__FILE__

        # do something only when this file is executed. 
        main(ARGS)
    end
end # module
