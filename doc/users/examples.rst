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
`Axes <http://matplotlib.sourceforge.net/api/axes_api.html>`__ instance method, 
with some extra pre/post processing and argument checking. 
You can also plot on the map directly with the matplotlib 
`pyplot <http://matplotlib.sourceforge.net/api/pyplot_api.html>`__ interface,
or the `OO api <http://matplotlib.sourceforge.net/examples/api/index.html>`__, 
using the `Axes <http://matplotlib.sourceforge.net/api/axes_api.html>`__ instance 
associated with the Basemap.

For more specifics of how to use the Basemap instance methods,
see :ref:`api-index`.

Here are the examples (many of which utilize the 
`netcdf4-python <http://netcdf4-python.googlecode.com>`__ module
to retrieve datasets over http):

* Plot contour lines on a basemap

.. literalinclude:: figures/contour1.py

.. image:: figures/contour1.png

* Plot precip with filled contours

.. literalinclude:: figures/plotprecip.py

.. image:: figures/plotprecip.png

* Plot sea-level pressure weather map with labelled highs and lows

.. literalinclude:: figures/plothighsandlows.py
 
.. image:: figures/plothighsandlows.png

* Plot hurricane tracks from a shapefile

.. literalinclude:: figures/hurrtracks.py
 
.. image:: figures/hurrtracks.png

* Plot etopo5 topography/bathymetry data as an image (with 
  and without shading from a specified light source).

.. literalinclude:: figures/plotetopo5.py
 
.. image:: figures/etopo5.png

.. image:: figures/etopo5_shaded.png

* Plot markers at locations of `ARGO <http://www.argo.ucsd.edu/>`__ floats.

.. literalinclude:: figures/plotargo.py
 
.. image:: figures/plotargo.png

* Pseudo-color plot of SST and sea ice analysis.

.. literalinclude:: figures/plotsst.py
 
.. image:: figures/plotsst.png

* Plotting wind vectors and wind barbs.

.. literalinclude:: figures/plotwindvec.py
 
.. image:: figures/plotwindvec1.png

.. image:: figures/plotwindvec2.png

* Draw great circle between NY and London.

.. literalinclude:: figures/plotgreatcircle.py
 
.. image:: figures/plotgreatcircle.png

* Draw day-night terminator on a map.

.. literalinclude:: figures/plotdaynight.py
 
.. image:: figures/plotdaynight.png
