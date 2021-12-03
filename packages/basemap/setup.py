#! /usr/bin/env python
# -*- coding: utf8 -*-
# flake8: noqa: E122
from __future__ import (absolute_import, division, print_function)

import io
import os
import sys
from setuptools import setup
from setuptools import find_packages
from setuptools.dist import Distribution
from setuptools.extension import Extension


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


# Define GEOS install directory (from environment variable or trying to guess).
geos_installdir = os.environ.get("GEOS_DIR", None)
if geos_installdir is not None:
    geos_version = checkversion(geos_installdir)
else:
    # Define some default locations to find GEOS.
    geos_search_locations = [
        os.path.expanduser("~/local"),
        os.path.expanduser("~"),
        "/usr/local",
        "/usr",
        "/opt/local",
        "/opt",
        "/sw"
    ]
    # Loop over the default locations to see if we find something.
    for folder in geos_search_locations:
        # sys.stdout.write('checking for GEOS lib in %s ....\n' % folder)
        geos_version = checkversion(folder)
        if geos_version is not None and geos_version >= '"3.1.1"':
            # sys.stdout.write(
            #     "GEOS lib (version %s) found in %s\n" %
            #     (geos_version[1:-1], folder))
            geos_installdir = folder
            break
    if geos_installdir is None:
        raise OSError("\n".join([
            "Cannot find GEOS library in standard locations ('{0}').",
            "Please install the corresponding packages using your",
            "software management system or set the environment variable",
            "GEOS_DIR to point to the location where GEOS is installed",
            "(for example, if 'geos_c.h' is in '/usr/local/include'",
            "and 'libgeos_c' is in '/usr/local/lib', then you need to",
            "set GEOS_DIR to '/usr/local'",
        ]).format("', '".join(geos_search_locations)))

# Define include dirs.
include_dirs = [
    os.path.join(geos_installdir, "include")
]
if not any("--" + opt in sys.argv for opt in Distribution.display_option_names +
           ["help-commands", "help"]) or sys.argv[1] == "egg_info":
    # Do not import numpy for just querying the package (taken from h5py setup).
    import numpy
    include_dirs.append(numpy.get_include())

# Define library dirs.
library_dirs = [
    os.path.join(geos_installdir, "lib"),
    os.path.join(geos_installdir, "lib64"),
]

# Define runtime library dirs.
# Don't use runtime_library_dirs on windows (workaround for a distutils bug):
#     http://bugs.python.org/issue2437)
if sys.platform == "win32":
    runtime_lib_dirs = []
else:
    runtime_lib_dirs = library_dirs

# Define `_geoslib` extension module. It cannot be installed in the
# `mpl_toolkits.basemap` namespace or `Basemap` objects will not be pickleable.
ext_modules = [
    Extension(**{
        "name":
            "_geoslib",
        "sources": [
            "src/_geoslib.c",
        ],
        "libraries": [
            "geos_c",
        ],
        "include_dirs":
            include_dirs,
        "library_dirs":
            library_dirs,
        "runtime_library_dirs":
            runtime_lib_dirs,
    }),
]

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
        ext_modules,
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
