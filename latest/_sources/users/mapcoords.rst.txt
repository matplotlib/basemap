.. _mapcoords:

Converting to and from map projection coordinates
=================================================

In order to plot data on a map, the coordinates of the data must
be given in map projection coordinates.
Calling a Basemap class instance with the arguments lon, lat will
convert lon/lat (in degrees) to x/y map projection coordinates
(in meters). The inverse transformation is done if the optional keyword
``inverse`` is set to True.  
Here's an example that uses this feature to plot a marker and some text to 
denote the location of Boulder, CO, given the lat/lon position.

.. plot:: users/figures/plotboulder.py

.. toctree::
