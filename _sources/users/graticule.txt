.. _graticule:

Drawing and Labelling Parallels and Meridians
=============================================

Most maps include a graticule grid, a reference network of labelled
latitude and longitude lines. Basemap does this with the 
:func:`~mpl_toolkits.basemap.Basemap.drawparallels` and
:func:`~mpl_toolkits.basemap.Basemap.drawmeridians` instance methods.
The longitude and latitude lines can be labelled where they intersect
the map projection boundary.  There are a few exceptions:  meridians
and parallels cannot be labelled on maps with 
``proj`` set to ``ortho`` (orthographic), ``geos`` (geostationary),
``vandg`` (van der Grinten) or ``nsper`` (near-sided perspective),
and meridians cannot be labelled on maps with 
``proj`` set to ``ortho`` (orthographic), ``geos`` (geostationary),
``vandg`` (van der Grinten), ``nsper`` (near-sided perspective),
``moll`` (Mollweide), ``hammer`` (Hammer), or ``sinu``
(sinusoidal).  This is because the lines can be very close 
together where they intersect the boundary on these maps, so that
they really need to be labelled manually on the interior of the plot.
Here's an example that shows how to draw parallels and meridians
and label them on different sides of the plot.

.. plot:: users/figures/graticule.py

.. toctree::
