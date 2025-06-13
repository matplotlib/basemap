.. _aeqd:

Azimuthal Equidistant Projection
================================

The shortest route from the center of the map
to any other point is a straight line in the azimuthal
equidistant projection. 
So, for the specified point, all points that lie on a circle around
this point are equidistant on the surface of the earth on this projection.
The specified point ``lon_0, lat_0`` shows up as a black dot in the center of the map.

Here's an example using the width and height keywords to specify the map region.

.. plot:: users/figures/aeqd.py

If both the width/height and corner lat/lon keywords are omitted, the whole world is 
plotted in a circle.

.. plot:: users/figures/aeqd_fulldisk.py
