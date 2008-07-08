.. _azeqd:

Azimuthal Equidistant Projection
================================

The shortest route from the center of the map
to any other point is a straight line in the azimuthal
equidistant projection. 
So, for the specified point, all points that lie on a circle around
this point are equidistant on the surface of the earth on this projection.
The specified point ``lon_0, lat_0`` shows up as a black dot in the center of the map.

.. literalinclude:: figures/azeqd.py

.. image:: figures/azeqd.png
