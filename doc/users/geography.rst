.. _geography:

Drawing a Map Background
========================

Basemap includes the GSSH (now
`GSHHG <https://www.soest.hawaii.edu/pwessel/gshhg/>`_)
coastline dataset, as well as datasets for rivers, state and
country boundaries from 
`GMT <http://gmt.soest.hawaii.edu>`_.
These datasets can be used to draw coastlines, rivers and political
boundaries on maps at several different resolutions.  The relevant Basemap 
methods are:

* :func:`~mpl_toolkits.basemap.Basemap.drawcoastlines`: draw coastlines.
* :func:`~mpl_toolkits.basemap.Basemap.fillcontinents`: color the interior
  of continents (by filling the coastline polygons).
  Unfortunately, the fillcontinents method doesn't always do the right thing.
  Matplotlib always tries to fill the inside of a polygon.  Under certain situations,
  what is the inside of a coastline polygon can be ambiguous, and the 
  outside may be filled instead of the inside.  
  In these situations, the recommended workaround is to use the 
  :func:`~mpl_toolkits.basemap.Basemap.drawlsmask` method to 
  overlay an image with different colors specified for land and water regions
  (see below).
* :func:`~mpl_toolkits.basemap.Basemap.drawcountries`: draw country boundaries.
* :func:`~mpl_toolkits.basemap.Basemap.drawstates`: draw state boundaries
  in North America.
* :func:`~mpl_toolkits.basemap.Basemap.drawrivers`: draw rivers.

Instead of drawing coastlines and political boundaries, an image can be
used as a map background.  Basemap provides several options for this:

* :func:`~mpl_toolkits.basemap.Basemap.drawlsmask`: draw a high-resolution 
  land-sea mask as an image, with land and ocean colors specified. The land-sea
  mask is derived from the GSHHS coastline data, and there are several 
  coastline options and pixel sizes to choose from.
* :func:`~mpl_toolkits.basemap.Basemap.bluemarble`: draw a NASA
  `Blue Marble <http://visibleearth.nasa.gov/view_set.php?categoryID=2363>`_
  image as a map background.
* :func:`~mpl_toolkits.basemap.Basemap.shadedrelief`: draw a  
  `shaded relief <http://www.shadedrelief.com>`_ image
  as a map background.
* :func:`~mpl_toolkits.basemap.Basemap.etopo`: draw an  
  `etopo <http://www.ngdc.noaa.gov/mgg/global/global.html>`_
  relief image as map background.
* :func:`~mpl_toolkits.basemap.Basemap.warpimage`: use an abitrary
  image as a map background.  The image must be global, covering the
  world in lat/lon coordinates from the international dateline eastward
  and the South Pole northward.

Here are examples of the various ways to draw a map background.

1. Draw coastlines, filling ocean and land areas.

.. plot:: users/figures/background1.py

2. Draw a land-sea mask as an image.

.. plot:: users/figures/background2.py

3. Draw the NASA 'Blue Marble' image.

.. plot:: users/figures/background3.py

4. Draw a shaded relief image.

.. plot:: users/figures/background4.py
 
5. Draw an etopo relief image.

.. plot:: users/figures/background5.py
