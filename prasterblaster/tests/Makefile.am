
AM_CLFLAGS =  -W -Wall -Wcast-align -Wstrict-prototypes -Wmissing \
		-prototypes -Wpointer-arith -Wshadow -Wsign-compare \
		-Wformat=2 -Wno-format-y2k -Wimplicit -Wmissing-braces \
		-Wnested-externs -Wparentheses -Wtrigraphs -g

noinst_PROGRAMS = tests systemtest
tests_CXXFLAGS  = -I../src/ -lgtest
tests_LDFLAGS = $(all_libraries) -lgdal -lgtest -lgtest_main -L../src/ -L../src/gtest/ -L../src/gdal/lib/
tests_CPPFLAGS = -I ../src/gtest/include -I../src/gdal/include
tests_LDADD = ../src/librasterblaster.la
tests_SOURCES = check-reprojection_tools.cc check-rastercoordtransformer.cc check_gdal.cc

systemtest_CXXFLAGS  = -I../src/ -lgtest
systemtest_LDFLAGS = $(all_libraries) -lgdal -lgtest -L../src/ -L../src/gtest/ -L../src/gdal/lib/
systemtest_CPPFLAGS = -I ../src/gtest/include -I../src/gdal/include
systemtest_LDADD = ../src/librasterblaster.la
systemtest_SOURCES = systemtest.cc ../src/demos/prasterblaster-pio.cc ../src/demos/sptw.cc

