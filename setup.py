from __future__ import (absolute_import, division, print_function)

import io
import os
import sys
from setuptools import setup
from setuptools import find_packages
from setuptools.dist import Distribution

if sys.version_info < (2, 6):
    raise SystemExit("""matplotlib and the basemap toolkit require Python 2.6 or later.""")

# Do not require numpy for just querying the package
# Taken from the netcdf-python setup file (which took it from h5py setup file).
inc_dirs = []
if not any('--' + opt in sys.argv for opt in Distribution.display_option_names +
           ['help-commands', 'help']) or sys.argv[1] == 'egg_info':
    import numpy
    # Use numpy versions if they are available.
    from numpy.distutils.core import setup, Extension
    # append numpy include dir.
    inc_dirs.append(numpy.get_include())


def get_content(name, splitlines=False):
    """Return the file contents with project root as root folder."""

    here = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(here, name)
    with io.open(path, encoding="utf-8") as fd:
        content = fd.read()
    if splitlines:
        content = [row for row in content.splitlines() if row]
    return content


def checkversion(GEOS_dir):
    """check geos C-API header file (geos_c.h)"""
    try:
        f = open(os.path.join(GEOS_dir, 'include', 'geos_c.h'))
    except IOError:
        return None
    geos_version = None
    for line in f:
        if line.startswith('#define GEOS_VERSION'):
            geos_version = line.split()[2]
    return geos_version

# get location of geos lib from environment variable if it is set.
if 'GEOS_DIR' in os.environ:
    GEOS_dir = os.environ.get('GEOS_DIR')
else:
# set GEOS_dir manually here if automatic detection fails.
    GEOS_dir = None

user_home = os.path.expanduser('~')
geos_search_locations = [user_home, os.path.join(user_home, 'local'),
                         '/usr', '/usr/local', '/sw', '/opt', '/opt/local']

if GEOS_dir is None:
    # if GEOS_dir not set, check a few standard locations.
    GEOS_dirs = geos_search_locations
    for direc in GEOS_dirs:
        geos_version = checkversion(direc)
        sys.stdout.write('checking for GEOS lib in %s ....\n' % direc)
        if geos_version is None or geos_version < '"3.1.1"':
            continue
        else:
            sys.stdout.write('GEOS lib (version %s) found in %s\n' %\
                    (geos_version[1:-1],direc))
            GEOS_dir = direc
            break
else:
    geos_version = checkversion(GEOS_dir)

if GEOS_dir is None:
    raise SystemExit("""
Can't find geos library in standard locations ('%s').
Please install the corresponding packages using your
systems software management system (e.g. for Debian Linux do:
'apt-get install libgeos-3.3.3 libgeos-c1 libgeos-dev' and/or
set the environment variable GEOS_DIR to point to the location
where geos is installed (for example, if geos_c.h
is in /usr/local/include, and libgeos_c is in /usr/local/lib,
set GEOS_DIR to /usr/local), or edit the setup.py script
manually and set the variable GEOS_dir (right after the line
that says "set GEOS_dir manually here".""" % "', '".join(geos_search_locations))
else:
    geos_include_dirs=[os.path.join(GEOS_dir,'include')] + inc_dirs
    geos_library_dirs=[os.path.join(GEOS_dir,'lib'),os.path.join(GEOS_dir,'lib64')]

# can't install _geoslib in mpl_toolkits.basemap namespace,
# or Basemap objects won't be pickleable.

# don't use runtime_library_dirs on windows (workaround
# for a distutils bug - http://bugs.python.org/issue2437).
if sys.platform == 'win32':
    runtime_lib_dirs = []
else:
    runtime_lib_dirs = geos_library_dirs

extensions = [ Extension("_geoslib",['src/_geoslib.c'],
                         library_dirs=geos_library_dirs,
                         runtime_library_dirs=runtime_lib_dirs,
                         include_dirs=geos_include_dirs,
                         libraries=['geos_c']) ]

# To create the source .tar.gz file:  python setup.py sdist
# To create the universal wheel file: python setup.py bdist_wheel --universal
setup(**{
    "name":
        "basemap",
    "version":
        "1.2.2+dev",
    "license":
        "OSI Approved",
    "description":
        "Plot data on map projections with matplotlib",
    "long_description": """
An add-on toolkit for matplotlib that lets you plot data on map
projections with coastlines, lakes, rivers and political boundaries.
See http://matplotlib.org/basemap/users/examples.html for examples of
what it can do.""",
    "long_description_content_type":
        "text/markdown",
    "url":
        "https://matplotlib.org/basemap",
    "author":
        "Jeff Whitaker",
    "author_email":
        "jeffrey.s.whitaker@noaa.gov",
    "maintainer":
        "Víctor Molina García",
    "maintainer_email":
        "molinav@users.noreply.github.com",
    "classifiers": [
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    "keywords": [
        "GIS",
        "maps",
        "plots",
    ],
    "namespace_packages": [
        "mpl_toolkits",
    ],
    "package_dir":
        {"": "lib"},
    "packages":
        find_packages(where="lib"),
    "ext_modules":
        extensions,
    "python_requires":
        ", ".join([
            ">=2.6",
            "!=3.0.*",
            "!=3.1.*",
            "<4",
        ]),
    "install_requires":
        get_content("requirements.txt", splitlines=True),
    "project_urls": {
        "Bug Tracker":
            "https://github.com/matplotlib/basemap/issues",
        "Documentation":
            "https://matplotlib.org/basemap/",
        "Source":
            "https://github.com/matplotlib/basemap",
    },
})
