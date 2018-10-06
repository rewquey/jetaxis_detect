dynlib
======

Dynlib comprises a collection of easy-to-use Fortran functions that provide 
interesting atmospheric and oceanographic diagnostics to be used with gridded
data sets. The collection is complemented by python utilities that further 
simplify the work with the Fortran functions.

Features
--------

- High-performance diagnostics and detection algorithms written in Fortran
- Input/Output library for netCDF, matlab and python data files
- A highly customisable plotting library based on matplotlib, that pads some of
  the sharp edges of matplotlib and provides a large collection of sensible 
  defaults for automatic plotting

Installation
------------

dynlib is available on all UiB Linux machines. Put the lines::
  export SHARED='/Data/gfi/users/local'
  export PATH="$PATH:$SHARED/bin"
  export PYTHONPATH="$SHAREDLIB/lib/python2.7/site-packages"
into your ``~/.bashrc`` to use it.

If you want to start developing new tools for dynlib, you will need to install
your personal version of the library. Follow the developer section in the 
dynlib documentation.

Contribute
----------

- Issue Tracker: github.com/dynlib/dynlib/issues
- Source Code: github.com/dynlib/dynlib

Support
-------

If you are having issues, please let us know: dynlib@gfi.uib.no.

License
-------

The project license is not determined yet.
