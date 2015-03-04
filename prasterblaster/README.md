
libRasterBlaster {#mainpage}
===============

    librasterblaster is made up of components that implement different
    parts of a parallel reprojection program. The components are
    designed so that they can be easily replaced with a user-provided
    component, hopefully allowing experimentation with new techniques.

    The sample implementation, demonstrating the use of
    librasterblaster is called \link prasterblasterpio \endlink.

Contents
--------
* libRasterBlaster setup and installation
* Background information and terminology
* librasterblaster API
* SPTW

libRasterBlaster Setup and Installation
---------------------------------------

### Pre-install steps

libRasterBlaster depends on GDAL so go to http://gdal.org and fetch
the sources or binaries for your system. Alternatively you can use
your system's package manager to install libgdal. Whatever method you
use ensure that libgdal is compiled with bigtiff and proj4 support.


### Compiling and Testing

Next check out the latest librasterblaster sources from the git
repository:

    git clone https://github.com/kornholi/prasterblaster.git

Then change to the prasterblaster directory.

    cd prasterblaster

To build minimal versions of GDAL and proj4 for use with
prasterblaster, use the buildgdal.sh script:

    bash buildgdal.sh

Run cmake to generate the makefile. The MPI implementation will be
automatically detected by cmake.

    mkdir build && cd build
    cmake ..

After running buildgdal.sh cmake should automatically
find the local gdal library.

After successfully running cmake you can build librasterblaster and
the demo program by calling make.

    make

To run the included tests run:

    make check

If any of the tests fail something is wrong and prasterblaster will
not function correctly.

To build the documentation use:

    make doc

The generated documentation can be found in the html/ directory.


libRasterBlaster now comes with some demo input files.

    bash ../tests/testdata/getdata.sh

The demo program is called 'prasterblasterpio'

    ./prasterblasterpio --t_srs +proj=moll ../tests/testdata/glc_geographic_30sec.tif ../ltests/testoutput/glc_mollweide_30sec.tif

You can also specify projections using EPSG

    ./prasterblasterpio --t_srs "EPSG:9842"

prasterblasterpio can be run in parallel as well:

    mpirun -n 100 ./prasterblasterpio --t_srs +proj=moll -n 4 tests/testdata/glc_geographic_30sec.tif tests/testoutput/glc_mollweide_30sec.tif


Background information and terminology
--------------------------------------

### Coordinate spaces

Most of the functionality in librasterblaster is concerned with the
manipulation of three distinct coordinate spaces that I call: raster,
projected, and geographic.

The raster space is used to address the individual pixel values in the
raster. It is a 2D cartesian coordinate plane with (0,0) in the
top-left corner. The y-axis grows down and the x-axis grows to the
right.

The projected coordinate space is a cartesian plane determined by the
map projection used. In projected coordinates (0,0) is in the center
of the space and the y-axis grows up.

Converting from projected coordinates to raster coordinates and back
is a simple transformation.

Geographic coordinates are the spheroidal longitude and latitude used
to describe locations on the earth. Converting projected coordinates
to geographic requires the use of map projection equations.



#### Projected

The projected coordinate system a raster uses is specified by its
projection string and its location in the projected space is specified
by its upper-left corner.

For example:

$ gdalinfo tests/testdata/nlcd2006_landcover_4-20-11_se5.tif |head -30

produces the output:


    Driver: GTiff/GeoTIFF
    Files: tests/testdata/nlcd2006_landcover_4-20-11_se5.tif
    Size is 161190, 104424
    Coordinate System is:
    PROJCS["Albers Conical Equal Area",
        GEOGCS["NAD83",
            DATUM["North_American_Datum_1983",
                SPHEROID["GRS 1980",6378137,298.2572221010002,
                    AUTHORITY["EPSG","7019"]],
                AUTHORITY["EPSG","6269"]],
            PRIMEM["Greenwich",0],
            UNIT["degree",0.0174532925199433],
            AUTHORITY["EPSG","4269"]],
        PROJECTION["Albers_Conic_Equal_Area"],
        PARAMETER["standard_parallel_1",29.5],
        PARAMETER["standard_parallel_2",45.5],
        PARAMETER["latitude_of_center",23],
        PARAMETER["longitude_of_center",-96],
        PARAMETER["false_easting",0],
        PARAMETER["false_northing",0],
        UNIT["metre",1,
            AUTHORITY["EPSG","9001"]]]
    Origin = (-2493045.000000000000000,3310005.000000000000000)
    Pixel Size = (30.000000000000000,-30.000000000000000)
    Metadata:
      AREA_OR_POINT=Area
    Image Structure Metadata:
      INTERLEAVE=BAND
    Corner Coordinates:
    Upper Left  (-2493045.000, 3310005.000) (130d13'58.18"W, 48d42'26.63"N)

The coordinate system is specified in the WKT format and the
upper-left corner of the raster is listed on the last line.

Using this information we could convert any raster coordinate to its
equivalent projection coordinate. For example:

#### Raster

Where is the raster coordinate (33, 150) in projected coordinates?

    X_proj = (X * pixel_size) + upper-left-x
           = (33 * 30) + -2493045
           = -2492055

    Y_proj = upper-left-y - (Y * pixel_size)
           = 3310005 - (150 * 30)
           = 326505

So the equivalent projected coordinate is (-2492055, 326505).


#### Geographic

Finding the geographic coordinate requires a more complicated transformation, so we'll use the invproj program:

    $ invproj +proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
    -249205  326505
    98d27'56.789"W  25d57'33.691"N <-- This is the geographic coordinate

libRasterBlaster API
--------------------

    librasterblaster is a library of tools to transform the map
    projection of raster datasets, potentially using many processors
    to speed up the computationally intensive processing required. The
    hope is that this will allow further research into I/O,
    partitioning, load balancing and other HPC concerns with regard to
    map reprojection.

    It is built on the MPI API for IPC and uses the GDAL library to
    access geospatial datasets. If you aren't familiar with the GDAL
    API check out this tutorial for a quick introduction:
    http://www.gdal.org/gdal_tutorial.html .

    The first is that you have a raster in a given projection and you
    want a new raster with a different map projection. The new output
    raster may be a very different shape than the input, so the extent
    of the output raster has to be calculated. This is done with the
    "ProjectedMinbox" function. It will return a rectangle, in the
    output projection's coordinate system, that will contain the input
    raster's area. This rectangular area may be proportioned very
    differently, depending on the projections involved.



    INPUT:
    +------------------------------------------------------------------+
    |                                       --+                        |
    |                                  \---/  |                        |
    |            /--------+            |\     |                        |
    |        +---         |            |      |    ---------\          |
    |        |            |             \      ---/          --        |
    |        \            |             +-                     -+      |
    |         \          /                \+                    \      |
    |         \          |                 |                     |     |
    |          X         |                 |                     |     |
    |         |          /                 |                     /     |
    |         |         /                   \                   |      |
    |        /          |                   |                   |      |
    |       /          /                    |                  /       |
    |      /         --                     |                //        |
    |     |       --/                       |            /---          |
    |     +\   --/                         /        /----              |
    |       --/                          \/   /-----                   |
    |                                     ----                         |
    +------------------------------------------------------------------+

                             +---+
                             |   |
                             |   |
                             |   |
                          -\-+---+--- 
                            -\   -/
                              --/

    OUTPUT:
    .......................-------------------........................
    .               ------/                   \------                .
    .           ---/  /  --            +-----\       \---            .
    .        --/     /     -\          \      \          \--         .
    .      -/       /        -\         |      \            \-       .
    .    -/         |          --       |     --\--------\    \-     .
    .  -/          /             \       \                \     \-   .
    . /           --              -\      -\              |       \  .
    . |             \--             \       \              |      |  .
    ./                 \+            X      |              |       \ .
    .|                  |          /-       |               \      | .
    .\                  |        /-        /+               /      / .
    . |                 |       /        /-                /      |  .
    . \                 /       /       /                --       /  .
    .  -\             -/       |       |              --/       /-   .
    .    -\         -/        /+       /           --/        /-     .
    .      -\     -/        /-        |         --/         /-       .
    .        --\----     /-/          |      --/         /--         .
    .           ---\\----              \  --/        /---            .
    .               ------\             \/    /------                .
    .......................-------------------........................



SPTW
----

The Simple Parallel Tiff Writer (SPTW) is a correct and almost optimal
implemention of parallel tiff file output using the MPI-IO api. The
MPI specification imposes three requirements to guarantee sequential
consistency:

    Case 1: $fh_1 FH_1$ All operations on fh1 are sequentially consistent
    if atomic mode is set. If nonatomic mode is set, then all operations
    on fh1 are sequentially consistent if they are either nonconcurrent,
    nonconflicting, or both.

    Case 2: $fh_1a FH_1$ and $fh_1b FH_1$ Assume A1 is a data access
    operation using fh1a, and A2 is a data access operation using fh1b. If
    for any access A1, there is no access A2 that conflicts with A1, then
    MPI guarantees sequential consistency.

    However, unlike POSIX semantics, the default MPI semantics for
    conflicting accesses do not guarantee sequential consistency. If A1
    and A2 conflict, sequential consistency can be guaranteed by either
    enabling atomic mode via the MPI_FILE_SET_ATOMICITY routine, or
    meeting the condition described in Case 3 below.

    Case 3: $fh_1 FH_1$ and $fh_2 FH_2$


SPTW provides sequential consistency only when all write operations
are nonconflicting, that is, each write operation accesses a distinct
section of the raster file. In the prasterblaterpio demo program each
pixel of the output raster is contained in only one partition and each
partition is only assigned to one process. This ensures all write
accesses are nonconflicting and that sequential consistency is
maintained.

librasterblaster algorithm
--------------------------

The librasterblaster uses a raster reprojection technique called
inverse reprojection. The output file is first created and then the
value of each output pixel is calculated.

### Serial Algorithm

#### 1. Calculate the size of the output raster in projected coordinates

First we iterate around the edges of the input raster, transforming
the coordinates of each pixel to the projection of the output
raster. A minbox is calculated as the smalled box in the projected
coordinates that would contain all of the calculated output points.

#### 2. Calculate the size of the output raster in raster coordinates

Using the size of the output pixels we then calculate the number of
rows and columns in the output raster like this:

columns = (lr_x - ul_x) / pixel_size
rows    = (ul_y - lr_y) / pixel_size

#### 3. Now create the output raster file

Create the output raster file with the user-provided projection and
the calculated raster size.

#### 4. Resample each output pixel

For each output pixel calculate the corresponding area in the input
raster and use the input raster values to calculate a new output pixel
value.

Write the calculated pixel value to the output file.

### Parallel Algorithm

The parallel algorithm is very similar to the serial algorithm except
now we resample the output raster in independent parts the can be
calculated in parallel.

#### 1. Calculate the size of the output raster in projected coordinates

First we iterate around the edges of the input raster, transforming
the coordinates of each pixel to the projection of the output
raster. A minbox is calculated as the smalled box in the projected
coordinates that would contain all of the calculated output points.

The librasterblaster function that performs this operation is called
librasterblaster::ProjectedMinbox.

#### 2. Calculate the size of the output raster in raster coordinates

Using the size of the output pixels we then calculate the number of
rows and columns in the output raster like this:

    columns = (lr_x - ul_x) / pixel_size
    rows    = (ul_y - lr_y) / pixel_size

#### 3. Now create the output raster file

Create the output raster file with the user-provided projection and
the calculated raster size.

Rasters can be created with the
librasterblaster::ProjectedRaster::CreateOutput function or with the
sptw::create_raster function.

#### 4. Partition the Output raster

Partition the output raster into n continuous parts. For each output
partition we then calculate the minbox of the corresponding input
area.

The librasterblaster::BlockPartition function can be used to partition
the output raster. Though many other partition techniques can be
applied.

To find the matching input raster area the function
librasterblaster::RasterMinbox applied with the output partition.


SPTW
----

SPTW is the Simple Parallel Tiff Writer. It's a library that uses the
MPI-IO api to allow multiple processes to consistently write a single
raster file. All of the processes in the process group must call the
sptw::open_raster function but once opened processes can make
independent calls to sptw::write_area. 

* As long as processes write to independent sections of the raster,
sequential consistency is preserved.

* All SPTW calls must be preceded by a call to MPI_Init().

* Every process in the MPI process group must open the same file at the same time.

API Examples
------------

Several examples of the API are provided. Check out
src/demos/prasterblaster-simple.cc for an example of a simple, serial
librasterblaster program or src/demos/prasterblaster-pio.cc for a
fully parallel implementation which uses SPTW.
