.. _examples:

Plotting data on a map (Example Gallery)
========================================

Following are a series of examples that illustrate how to use
Basemap instance methods to plot your data on a map.  More examples
are included in the examples directory of the basemap source distribution.
There are a number of Basemap instance methods for plotting data:

* :func:`~mpl_toolkits.basemap.Basemap.contour`: draw contour lines.
* :func:`~mpl_toolkits.basemap.Basemap.contourf`: draw filled contours.
* :func:`~mpl_toolkits.basemap.Basemap.imshow`: draw an image.
* :func:`~mpl_toolkits.basemap.Basemap.pcolor`: draw a pseudocolor plot.
* :func:`~mpl_toolkits.basemap.Basemap.pcolormesh`: draw a pseudocolor plot (faster version for regular meshes).
* :func:`~mpl_toolkits.basemap.Basemap.plot`: draw lines and/or markers.
* :func:`~mpl_toolkits.basemap.Basemap.scatter`: draw points with markers.
* :func:`~mpl_toolkits.basemap.Basemap.quiver`: draw vectors.
* :func:`~mpl_toolkits.basemap.Basemap.barbs`: draw `wind barbs <http://en.wikipedia.org/wiki/Station_model#Plotted_winds>`__.
* :func:`~mpl_toolkits.basemap.Basemap.drawgreatcircle`: draw a `great circle <http://en.wikipedia.org/wiki/Great_circle>`__.

Many of these instances methods simply forward to the corresponding matplotlib
`Axes <https://matplotlib.org/api/axes_api.html>`__ instance method, 
with some extra pre/post processing and argument checking. 
You can also plot on the map directly with the matplotlib 
`pyplot <https://matplotlib.org/api/pyplot_api.html>`__ interface,
or the `OO api <https://matplotlib.org/examples/api/index.html>`__, 
using the `Axes <https://matplotlib.org/api/axes_api.html>`__ instance 
associated with the Basemap.

For more specifics of how to use the Basemap instance methods,
see :ref:`api-index`.

Here are the examples (many of which utilize the 
`netcdf4-python <http://unidata.github.io/netcdf4-python/>`__ module
to retrieve datasets over http):

* Plot contour lines on a basemap

.. plot:: users/figures/contour1.py

* Plot precip with filled contours

.. plot:: users/figures/plotprecip.py

* Plot sea-level pressure weather map with labelled highs and lows

.. plot:: users/figures/plothighsandlows.py
 
* Plot hurricane tracks from a shapefile

.. plot:: users/figures/hurrtracks.py
 
* Plot etopo5 topography/bathymetry data as an image (with 
  and without shading from a specified light source).

.. plot:: users/figures/plotetopo5.py
 
* Plot markers at locations of `ARGO <http://www.argo.ucsd.edu/>`__ floats.

.. plot:: users/figures/plotargo.py
 
* Pseudo-color plot of SST and sea ice analysis.

.. plot:: users/figures/plotsst.py
 
* Plotting wind vectors and wind barbs.

.. plot:: users/figures/plotwindvec.py
 
* Draw great circle between NY and London.

.. plot:: users/figures/plotgreatcircle.py
 
* Draw day-night terminator on a map.

.. plot:: users/figures/plotdaynight.py
