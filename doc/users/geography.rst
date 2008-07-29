.. _geography:

Drawing Coastlines, Rivers and Political Boundaries
===================================================

Basemap includes the 
`GSSH <http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html>`_
coastline dataset, as well as datasets for rivers, state and
country boundaries from 
`GMT <http://gmt.soest.hawaii.edu>`_.
These datasets can be used to draw coastlines, rivers and political
boundaries on maps at several different resolutions.  The relevant Basemap 
methods are:

* :func:`~mpl_toolkits.basemap.Basemap.drawcoastlines`: draw coastlines.
* :func:`~mpl_toolkits.basemap.Basemap.fillcontinents`: color the interior
  of continents (by filling the coastline polygons).
* :func:`~mpl_toolkits.basemap.Basemap.drawcountries`: draw country boundaries.
* :func:`~mpl_toolkits.basemap.Basemap.drawstates`: draw state boundaries
  in North America.
* :func:`~mpl_toolkits.basemap.Basemap.drawrivers`: draw rivers.

.. toctree::
