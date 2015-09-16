=================================================
Python Spatial Analysis Library
=================================================

.. Contents::

What is PySAL
--------------

PySAL is an open source cross-platform library of spatial analysis functions
written in Python. It is intended to support the development of high level
applications for spatial analysis.

It is important to underscore what PySAL is, and is not, designed to do. First
and foremost, PySAL is a library in the fullest sense of the word. Developers
looking for a suite of spatial analytical methods that they can incorporate
into application development should feel at home using PySAL. Spatial analysts
who may be carrying out research projects requiring customized scripting,
extensive simulation analysis, or those seeking to advance the state of the art
in spatial analysis should also find PySAL to be a useful foundation for their
work.

End users looking for a user friendly graphical user interface for spatial
analysis should not turn to PySAL directly. Instead, we would direct them to
projects like STARS and the GeoDaX suite of software products which wrap PySAL
functionality in GUIs. At the same time, we expect that with developments such
as the Python based plug-in architectures for QGIS, GRASS, and the toolbox
extensions for ArcGIS, that end user access to PySAL functionality will be
widening in the near future.

What is pPySAL
--------------
pPySAL is an extension to PySAL where select analytical methods are being
developed for parallel (multi-core) computation. pPySAL requires an 
installation of the core PySAL library in addition to all of the PySAL 
dependencies. 

pPySAL package structure
-----------------------

Currently PySAL consists of the following files and directories:

  LICENSE.txt
    PySAL license.

   fj_refactored.py
    A parallel Fisher-Jenks optimal map classifier.
   fj_test.py
    A test script to run the parllel Fisher-Jenks algorithm over
     a range of sample sizes using a variable number of processing
     cores

Website
-------
All things PySAL can be found here
    http://pysal.org/

Mailing Lists
-------------
Please see the developer's list here
    http://groups.google.com/group/pysal-dev

Help for users is here
    http://groups.google.com/group/openspace-list

Bug reports
-----------
To search for or report bugs, please see
    http://github.com/pysal/pPysal/issues

License information
-------------------
See the file "LICENSE.txt" for information on the history of this
software, terms & conditions for usage, and a DISCLAIMER OF ALL
WARRANTIES.
