module jtsk2etrs    

#
# cartographic trasnformation
#

include("coordinate_transformation/coordinate_transformation.jl")
include("coordinate_transformation/cartographic_projection.jl")
include("coordinate_transformation/ellipsoid.jl")
include("coordinate_transformation/etrs2jtsk.jl")

#
# data handling 
#

include("data_handling/geoid_model.jl")
include("data_handling/load_data.jl")
include("data_handling/table.jl")


if abspath(PROGRAM_FILE) == @__FILE__

    # do something only this file is executed. 
    do_something()

end

end # module
