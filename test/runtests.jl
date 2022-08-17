using Test
using jtsk2etrs

# Simple macro to allow named tests
macro namedtest(name, test)
    esc(:(@testset $name begin @test $test end))
end

@testset "Coordinate transformations" begin
    ###### Testing the relation between xyz and blh ######
    # Single point solutions 
    b, l, h = xyz2blh(4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6, wgs84)
    x, y, z = blh2xyz(49.20589174, 16.59283443, 324.272, wgs84)

    @namedtest "xyz2blh::scalar on wgs84 ellipsoid" isapprox([b,l,h], [49.20589174, 16.59283443, 324.272], atol=0.0000001)
    @namedtest "blh2xyz::scalar on wgs84 ellipsoid" isapprox([x,y,z], [4.0014706009668424e6, 1.192345306410758e6, 4.805795321625611e6], atol=1e-9)
        
    ### Vector solutions
    # input values
    bv::Vector{Float64} = [41.12345678; 67.321654987; -41.12345678; -67.321654987]
    lv::Vector{Float64} = [15.32145679;179.235468971; 245.79842613; 347.321654987]
    hv::Vector{Float64} = [321.4596; -123.5469; 8212.6544; -2354.1232]
    # expected values
    xres::Vector{Float64} = [4.640807730355608e6, -2.465908053190801e6, -1.9750362269111197e6, 2.4051583575709253e6]
    yres::Vector{Float64} = [1.2714495328089027e6, 32906.003486815716, -4.394332066337075e6, -541070.7280088026]
    zres::Vector{Float64} = [4.1729723702383414e6, 5.862223001747927e6, -4.178162280416622e6, -5.860164884959347e6]

    xv,yv,zv = blh2xyz(bv,lv,hv, wgs84)
    bbv,llb,hhv = xyz2blh(xres,yres,zres, wgs84)

    @namedtest "xyz2blh::vector on wgs84 ellipsoid" isapprox([xv,yv,zv], [xres,yres,zres], atol=1e-9)
    @namedtest "blh2xyz::vector on wgs84 ellipsoid" isapprox([bbv,llb,hhv], [bv,lv,hv], atol=1e-8)

    ###### Testing the relation between xyz and ruv ######
    ### Single point solutions     
    r, u, v = xyz2sph(x, y, z)

    rres, ures, vres = (6.3662487363816e6, 49.015454682101705, 16.59283443)
    xsph, ysph, zsph = sph2xyz(r,u,v)

    @namedtest "xyz2sph::scalar" isapprox([r, u, v], [rres, ures, vres], atol=1e-9)
    @namedtest "sph2xyz::scalar" isapprox([x,y,z], [xsph, ysph, zsph], atol=1e-9)
    ### Vector solutions

    rv, uv, vv  = xyz2sph(xres, yres, zres)
    xvsph, yvsph, zvsph = sph2xyz(rv, uv, vv)

    @namedtest "sph2xyz::vector" isapprox([xvsph, yvsph, zvsph], [xres,yres,zres], atol=1e-8)

end

@testset "Testing the cartographic projection for JTSK and JTSK05" begin
    b, l = jtsk2bl(712201.2924148451, 1.5644015536543208e6)

    # @namedtest "Converting JTSK yx to geodetic bl" isapprox([b,l], [45.15, 15.76823], atol=1e-5)
end

@testset "3D transformation between the WGS84 ellipsoid and Bessel's" begin
    btubo::Float64 = 49.2058916
    ltubo::Float64 = 16.592834502777777    
    htubo::Float64 = 324.374

    xtubo_wgs, ytubo_wgs, ztubo_wgs = blh2xyz(btubo, ltubo, htubo, wgs84)
    xtubo_bes, ytubo_bes, ztubo_bes = etrf2bessel(xtubo_wgs, ytubo_wgs, ztubo_wgs)
    xtubo_etrf, ytubo_etrf, ztubo_etrf = bessel2etrf(xtubo_bes, ytubo_bes, ztubo_bes)

    @namedtest "Tubo 'xyz' coordinates etrf2bessel::scalar" isapprox([xtubo_bes, ytubo_bes, ztubo_bes],[4.0020623292969726e6, 1.1924208321677265e6, 4.8062734167892095e6], atol=1e-5 )
    @namedtest "Tubo 'xyz' reverse transformation bessel2etrf::scalar" isapprox([xtubo_etrf, ytubo_etrf, ztubo_etrf ],[xtubo_wgs, ytubo_wgs, ztubo_wgs], atol=1e-3 )

    # Czepos coordinates
    bvec::Vector{Float64} = [48.96763098; 49.44590371; 49.68479282; 48.84962013; 50.23282443; 50.54027399; 50.41287907; 50.03954777; 49.01430046; 50.12522950; 49.68779800; 50.10238943; 49.96476928; 49.75782059; 49.40979677; 50.56264359; 49.33801148; 48.86403504; 49.91370184; 49.72657611; 50.35015140; 49.20589174; 49.83351384; 49.14800876; 49.39355391; 49.29748084; 50.23259320; 50.77170209]
    lvec::Vector{Float64} = [14.47527113; 12.92410465; 18.35318282; 17.12906703; 17.20816619; 14.14025005; 14.90596422; 15.78324593; 13.99593811; 14.45605687; 13.99825781; 13.72923698; 16.98095702; 16.47130791; 14.68021872; 15.9084478;  17.99101296; 16.03917764; 14.78561823; 13.35183533; 16.32224338; 16.59283443; 18.16383295; 15.00880698; 15.58639573; 17.40014282; 12.84187541; 15.05989126]
    hvec::Vector{Float64} = [456.223; 519.603;373.59 ;228.387;495.223;243.275;303.462;283.27 ;645.39 ;356.025;583.676;381.87 ;378.365;498.442;496.233;478.595 ;407.325;373.84 ;592.602;425.23 ;791.704;324.272;340.895;543.522;576.839;258.576;446.082;448.35 ]

    # conversion blh2xyz on wgs84
    xvec, yvec, zvec = blh2xyz(bvec, lvec, hvec, wgs84)

    xbes, ybes, zbes = etrf2bessel(xvec, yvec, zvec)
    xetrf, yetrf, zetrf = bessel2etrf(xbes, ybes, zbes)
    xbes2, ybes2, zbes2 = etrf2bessel(xetrf, yetrf, zetrf )

    @namedtest "CZEPOS::vector coordinates etrf   >> bessel >> etrf"   isapprox([xvec, yvec, zvec], [xetrf, yetrf, zetrf], atol=1e-3)
    @namedtest "CZEPOS::vector coordinates bessel >> etrf   >> bessel" isapprox([xbes, ybes, zbes], [xbes2, ybes2, zbes2], atol=1e-3)
end

@testset "Geodetic computations on the sphere" begin
    b1::Float64 = 49.2058916 # Tubo coordinates
    l1::Float64 = 16.592834502777777     # Tubo coordinates
    h1::Float64 = 324.374 # Tubo coordinates

    b2::Float64 = 39.208916 # 
    l2::Float64 = 10.592646 # 
    h2::Float64 = 212.32100 # 

    b3  ::Float64 = 44.208916 # 
    l3  ::Float64 = 15.511246 # 
    h3  ::Float64 = -112.3210 # 

end