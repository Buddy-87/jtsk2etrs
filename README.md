# The coordinate conversion between the ETRF2000 and S-JTSK-05

This project is based on [Metodika převodu mezi ETRF2000 a S-JTSK](https://www.cuzk.cz/Zememerictvi/Geodeticke-zaklady-na-uzemi-CR/GNSS/Nova-realizace-systemu-ETRS89-v-CR/Metodika-prevodu-ETRF2000-vs-S-JTSK-var2(101208).aspx). The methodology describes in detail the conversion between ETRF2000 coordinate systems and S-JTSK. (Note: According to the EUREF-TWG recommendation, it is intended for the implementation of the European reference frame to use the designation ETRF2000 instead of the previous ETRF2000(R05)). The intermediate variable for this conversion is the S-JTSK/05 system. The advent of Global Navigation Satellite Systems (GNSS) technology in the 1980s. century marked a revolutionary leap in the speed and, to a large extent, in the accuracy of positioning with "geodetic accuracy" (units of cm). In addition to the requirement to integrate geodetic foundations into of the pan-European coordinate system, it was necessary to enable the spatial connection
coordinate system with which GNSS technology works, with a plane system coordinates in the cartographic representation with which classic geodetic methods work. The implemented S-JTSK/05 coordinate system meets the following requirements:
- the system is primarily implemented by 3141 points of the DOPNUL network and "Selective maintenance"
- uses the spatial geodetic reference system for measurements by GNSS technology
framework in the implementation of 2005 – i.e. ETRF2000 – thereby guaranteeing continuity with
European spatial coordinate system
- uses a modified coordinate system for measurements in plane coordinates
Křovák's display - this guarantees continuity with the existing binding
S-JTSK coordinate system
- for the transformation of ellipsoidal heights (relative to the GRS80 ellipsoid) to heights
above sea level in the binding "Balt after alignment" system, a quasi-geoid model is used
CR-2005
- between coordinates in the ETRF2000 coordinate frame and plane coordinates
in the modified Křovák representation, an exact mathematical relationship applies
- the root mean square value of the position deviation between the coordinates in S-JTSK and SJTSK/05 is 13.3 cm, for work requiring less than 0.5 m accuracy, both systems are
pronoun


**Note: The methodology contains some mistakes (such as swaped values for 3D transformation between the Bessel's ellipsoid and the ellipsoid GRS80. The source code should contain correct math expressions etc.**
## How to install

To clone project:
```
git clone https://github.com/Buddy-87/jtsk2etrs.jl.git
```

## Run unit tests


You can run the unit tests to be sure that the package runs and computes stuff correctly.
It runs all the tests in folder ```./test/runtests.jl```

```bash
# navigate to ./test/..
julia runtests.jl 
```


## Run examples (CLI)

It is possible to use the Command line interface (CLI) for quick computation.

- `--file,-f`   (**mandatory**) - file containing the coordinates fo the points. See files *czepos_jtsk.txt*, *czepos_etrs.txt* as examples.
- `--transformation", -t`   (**mandatory**) - the direction of the transformation. Either from ETRF frame to JTSK or the opposite.
- `--output, -o`    (**default value is `result.txt`**) - name of the output.
- `--table` (**default value is `v1202`**) - table with the additional corrections for **y** and **x** coordinates in the JTSK reference frame.

```bash
julia jtsk2etrs.jl -f "../data/czepos_jtsk.txt" -t jtsk2etrf --table v1710
```

```bash
julia jtsk2etrs.jl -f "../data/czepos_etrs.txt" -t etrf2jtsk --table v1005
```

## References
- [Metodika převodu mezi ETRF2000 a S-JTSK](https://www.cuzk.cz/Zememerictvi/Geodeticke-zaklady-na-uzemi-CR/GNSS/Nova-realizace-systemu-ETRS89-v-CR/Metodika-prevodu-ETRF2000-vs-S-JTSK-var2(101208).aspx). 
- [The Julia Programming Language](https://julialang.org/)
- [Program for coordinate transformation between the coordinate systems used in the Czech Republic ](http://gisak.vsb.cz/GISacek/GISacek_2001/sbornik/Hanzlova/Hanzlova.htm)

