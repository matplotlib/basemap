Installation
============

Installing from PyPI
--------------------

Precompiled binary wheels for Windows and GNU/Linux are available in
PyPI (architectures x86 and x64, Python 2.7 and 3.5+) and can be
installed with `pip`_:

.. code-block:: sh

   python -m pip install basemap

Installing ``basemap`` will also install ``basemap-data``, containing the
minimal data assets required by ``basemap``. If you also need the
high-resolution data assets, you can install them with `pip`_ too:

.. code-block:: sh

   python -m pip install basemap-data-hires

Installing from conda-forge
---------------------------

For Miniforge users, ``basemap`` packages are available through the
``conda-forge`` channel for Windows and GNU/Linux (x64) as well as
for MacOS (x64 and arm64):

.. code-block:: sh

    conda install -c conda-forge basemap

Similarly to the PyPI installation, the high-resolution data assets
can be installed separately if needed:

.. code-block:: sh

    conda install -c conda-forge basemap-data-hires

Installing from source
----------------------

Optionally, you can also install ``basemap`` from its source hosted
on GitHub as indicated in the following steps:

1. Install pre-requisite Python modules:

   - `cython`_
   - `numpy`_

2. Download the ``basemap`` source code and move to the
   ``packages/basemap`` folder:

   .. code-block:: sh

      git clone --depth 1 https://github.com/matplotlib/basemap.git
      cd basemap/packages/basemap

3. Build the `GEOS`_ library. You may use the helper provided in the
   ``utils`` folder (please note that you need `CMake`_ and a working
   C compiler in advance):

   .. code-block:: sh

      export GEOS_DIR=<your desired location>
      python -c "import utils; utils.GeosLibrary('3.6.5').build(installdir='${GEOS_DIR}')"

   or you can link directly to the system library if it is already
   installed. ``GEOS_DIR`` must point to the GEOS installation prefix;
   e.g. if ``libgeos_c.so`` is located in ``/usr/lib`` and ``geos_c.h``
   is located in ``/usr/include``, then you must set ``GEOS_DIR`` to
   ``/usr``.

4. Build and install the ``basemap`` binary wheel:

   .. code-block:: sh

      python -m pip install .

   On GNU/Linux, if your Python was installed through a package
   management system, make sure that you have the Python header
   ``Python.h`` required to build Cython extensions (e.g. on
   Debian-like systems, you should have the package ``python-dev``
   installed).

5. Check that the package was installed correctly by executing:

   .. code-block:: sh

      python -c "from mpl_toolkits.basemap import Basemap"


.. _pip: https://pip.pypa.io/
.. _cython: https://github.com/cython/cython
.. _numpy: https://github.com/numpy/numpy
.. _GEOS: https://github.com/libgeos/geos
.. _CMake: https://cmake.org/
