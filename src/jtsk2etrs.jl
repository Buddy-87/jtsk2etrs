module jtsk2etrs    
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
    export Ellipsoid, grs80, wgs84
    export get_transformation_parameters, bessel2etrf, etrf2bessel


    #
    # data handling 
    #

    include("data_handling/geoid_model.jl")
    include("data_handling/load_data.jl")
    include("data_handling/table.jl")

    export Geoid, interpolate
    export load_table
    export Table, check_dimensions

    #
    # data plot
    #
    include("plot/plot_raw_data.jl")

    function main()::Integer
        println("Executing the program \"jtsk2etrs\".")

        return PROGRAM_STATUS_OK
    end


    if abspath(PROGRAM_FILE) == @__FILE__

        # do something only this file is executed. 
        main()
    end

end # module
