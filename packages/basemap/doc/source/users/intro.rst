Introduction
============

The matplotlib basemap toolkit is a library for plotting 2D data on maps
in `Python <http://www.python.org>`_. It is similar in functionality to
the `matlab mapping toolbox <http://www.mathworks.com/access/helpdesk/help/toolbox/map/map.shtml>`_,
the `IDL mapping facilities <http://www.msi.umn.edu/software/idl/tutorial/idl-mapping.html>`_, 
`GrADS <http://www.iges.org/grads/downloads.html>`_, or the 
`Generic Mapping Tools <http://gmt.soest.hawaii.edu/>`_. 
`PyNGL <http://www.pyngl.ucar.edu/>`_ and
`CDAT <http://www-pcmdi.llnl.gov/software/cdat/support/vcs/vcs.html>`_
are other libraries that provide similar capabilities in Python.

Basemap does not do any plotting on it's own, but provides the facilities to transform coordinates to one of 25 different map projections (using the 
`PROJ.4 <http://trac.osgeo.org/proj/>`_ C library).  `Matplotlib
<https://matplotlib.org>`_ is then
used to plot contours, images, vectors, lines or points
in the transformed coordinates.
Shoreline, river and political boundary
datasets (from `Generic Mapping Tools <http://gmt.soest.hawaii.edu/>`_)
are provided, along with methods for plotting them. The `GEOS library 
<http://geos.refractions.net>`_ is used internally to clip the coastline and polticial boundary features to the desired map projection region.

Basemap is geared toward the needs of earth scientists, particularly
oceanographers and meteorologists.  Jeff Whitaker originally wrote Basemap
to help in his research (climate and weather forecasting), since at the time 
`CDAT <http://www-pcmdi.llnl.gov/software/cdat/support/vcs/vcs.html>`_ was 
the only other tool in python for plotting data on map projections.  Over
the years, the capabilities of Basemap have evolved as scientists in other
disciplines (such as biology, geology and geophysics) requested and 
contributed new features.
