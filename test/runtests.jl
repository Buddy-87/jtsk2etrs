using Test

@testset "Coordinate transformation" begin
    @test begin
        x,y,z = xyz2blh(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6, wgs84)
        (49.20589174, 16.59283443, 324.272) 
        x == 49.20589174
        
    end
end
@test true