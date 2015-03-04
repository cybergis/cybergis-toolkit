#CyberGIS Toolkit

##Introduciton

CyberGIS Toolkit is a suite of loosely coupled open-source geospatial software components that provide computationally scalable spatial analysis and modeling capabilities enabled by advanced cyberinfrastructure. CyberGIS Toolkit represents a deep approach to CyberGIS software integration research and development and is one of the three key pillars of the CyberGIS software environment, along with [CyberGIS Gateway](http://gateway.cigi.illinois.edu/home/) and [GISolve Middleware](http://www.cigi.illinois.edu/dokuwiki/doku.php/projects/gisolve/index). The integration approach to building CyberGIS Toolkit is focused on developing and leveraging innovative computational strategies needed to solve computing- and data-intensive geospatial problems by exploiting high-end cyberinfrastructure resources such as supercomputing resources provided by the NSF Extreme Science and Engineering Discovery Environment ([XSEDE](http://xsede.org/)) and high-throughput computing resources on the Open Science Grid ([OSG](http://opensciencegrid.org/)). 

A rigorous process of software engineering and computational intensity analysis is applied to integrate an identified software component into the toolkit, including software building, testing, packaging, scalability and performance analysis, and deployment. This process includes three major steps:
Local build and test by software researchers and developers using continuous integration software or services such as [Travis CI](http://travis-ci.org/);
Continuous integration testing, portability testing, small-scale scalability testing on the National Middleware Initiative ([NMI](http://batlab.org/)) build and test facility; and
XSEDE-based evaluation and testing of software performance, scalability, and portability. By leveraging the high-performance computing expertise in the integration team of the NSF CyberGIS Project, large-scale problem-solving tests are conducted on various supercomputing environments on XSEDE to identify potential computational bottlenecks and achieve maximum problem-solving capabilities of each software installation.

##Software Components

1. [TauDEM](https://github.com/dtarb/TauDEM). [TauDEM](http://hydrology.usu.edu/taudem/) (Terrain Analysis Using Digital Elevation Models) is a suite of high-performance Digital Elevation Model (DEM) tools for watershed delineation and the extraction and analysis of hydrologic information from topography as represented by a Digital Elevation Model.
2. [Parallel PySAL](http://cybergis.cigi.uiuc.edu/cyberGISwiki/doku.php/ct/ppysal). The Parallel [PySAL](http://pysal.org/) library provides a set of scalable PySAL functions. Currently, parallel PySAL components implemented using the multiprocessing python library. The Fisher-Jenks classification algorithm is integrated in the toolkit;
3. [pRasterBlaster](http://cybergis.cigi.uiuc.edu/cyberGISwiki/doku.php/ct/prasterblaster). pRasterBlaster is a high-performance map reprojection software contributed by the Center of Excellence for Geospatial Information Science ([CEGIS](http://cegis.usgs.gov/)) within the [U.S. Geological Survey](http://usgs.gov/);
4. [PGAP](http://cybergis.cigi.uiuc.edu/cyberGISwiki/doku.php/ct/pgap). PGAP is a scalable Parallel Genetical Algorithm (PGA) solver for the Generalized Assignment Problem (GAP). This code provides an efficient PGA implementation for combinatorial optimization problem-solving and scaled up to 262K processor cores on BlueWaters with marginal communcation cost.
5. [Parallel Agent-Based Modeling (PABM)](http://cybergis.cigi.uiuc.edu/cyberGISwiki/doku.php/ct/pabm). PABM is an illustrative software for scalable spatially explicit agent-based modeling (ABM);
##Software Components Being Integrated

ScalaGAHealth (Scalable Geographic Analytics for Health-Related Data)
SpatialEvo (A high-performance evolutionary computing library for solving complex spatial optimization problems)

##Copyright and License

Each software component integrated in CyberGIS Toolkit is open source and has its own copyright and license.
- TauDEM: Please refer to [TauDEM](http://hydrology.usu.edu/taudem/) website
- Parallel PySAL: [BSD](http://opensource.org/licenses/BSD-3-Clause)
- pRasterBlaster: [GPL v3](http://opensource.org/licenses/gpl-3.0)
- PGAP: [NCSA open source license](http://opensource.org/licenses/NCSA)
- PABM: [NCSA open source license](http://opensource.org/licenses/NCSA)

##Contact

If you have any questions about CyberGIS Toolkit, please contact CyberGIS Helpdesk (help@cybergis.org).
