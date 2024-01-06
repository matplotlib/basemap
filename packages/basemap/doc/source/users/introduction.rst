Introduction
============

The matplotlib basemap toolkit is a library for plotting 2D data on maps
in `Python`_. It is similar in functionality to `GrADS`_, `GMT`_, the
`MATLAB Mapping Toolbox`_ and the `IDL Mapping Facilities`_. `CDAT`_
and `PyNGL`_ are other Python libraries with similar capabilities.

Basemap does not plot on its own, but provides the facilities to
transform coordinates to one of 25 different map projections (using
`pyproj`_ and therefore the `PROJ`_ C library). Then `matplotlib`_ is
used to plot contours, images, vectors, lines or points in the
transformed coordinates. Shoreline, river and political boundary
datasets (extracted from `GMT`_) are provided, together with methods
for plotting them. The `GEOS`_ library is used internally to clip the
coastline and political boundary features to the map projection region.

Basemap is geared towards the needs of Earth scientists, particularly
oceanographers and meteorologists. Jeff Whitaker originally wrote
Basemap to help in his research (climate and weather forecasting),
since at the time `CDAT`_ was the only other tool in Python for
plotting data on map projections. Over the years, the capabilities
of basemap have evolved as scientists in other disciplines (such as
biology, geology and geophysics) requested and contributed new features.


.. _Python: https://www.python.org/
.. _GMT: https://www.generic-mapping-tools.org/
.. _GrADS: http://cola.gmu.edu/grads/
.. _MATLAB Mapping Toolbox: https://www.mathworks.com/help/map/map.html
.. _IDL Mapping Facilities: https://www.nv5geospatialsoftware.com/docs/mapping_funct_list.html
.. _CDAT: https://cdat.llnl.gov/
.. _PyNGL: https://www.pyngl.ucar.edu/

.. _pyproj: https://pyproj4.github.io/pyproj
.. _PROJ: https://proj.org/
.. _matplotlib: https://matplotlib.org/
.. _GEOS: https://libgeos.org/
