.. _mapsetup:

Setting up the map
==================

In order to represent the curved surface of the earth on a two-dimensional
map, a map projection is needed. Since this cannot be done without
distortion, there are many map projections, each with it's own advantages
and disadvantages. Basemap provides 24 different map projections.
Some are global, some can only represent a portion of the globe. When
a :class:`~mpl_toolkits.basemap.Basemap` class instance is
created, the desired map projection must
be specified, along with information about the portion of the earth's
surface that the map projection will describe. There are two basic
ways of doing this. One is to provide the latitude and longitude values
of each of the four corners of the rectangular map projection region.
The other is to provide the lat/lon value of the center of the map
projection region along with the width and height of the region in
map projection coordinates. 

The class variable ``supported_projections`` is a dictionary containing 
information about all the projections supported by Basemap.  The keys
are the short names (used with the ``projection`` keyword to define
a projection when creating a ``Basemap`` class instance), and the values
are longer, more descriptive names.  The class variable ``projection_params``
is a dictionary that provides a list of parameters that can be used to
define the properties of each projection.  Following are examples that 
illustrate how to set up each of the supported projections. Note that
many map projection possess one of two desirable properties - they can be
equal-area (the area of features is preserved) or conformal (the shape of
features is preserved).  Since no map projection can have both at the same
time, many compromise between the two.

.. toctree::

    aeqd.rst
    gnomon.rst
    ortho.rst
    geos.rst
    nsper.rst
    moll.rst
    hammer.rst
    robin.rst
    eck4.rst
    kav7.rst
    mbtfpq.rst
    sinu.rst
    cyl.rst
    cass.rst
    merc.rst
    tmerc.rst
    omerc.rst
    poly.rst
    mill.rst
    gall.rst
    cea.rst
    lcc.rst
    laea.rst
    stere.rst
    eqdc.rst
    aea.rst
    pstere.rst
    plaea.rst
    paeqd.rst
    vandg.rst
