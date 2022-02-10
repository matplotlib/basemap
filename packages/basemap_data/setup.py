#! /usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa: E122
"""basemap_data -- Data assets for matplotlib basemap."""

import io
import os
import itertools
from setuptools import setup
from setuptools import find_packages


def get_content(name, splitlines=False):
    """Return the file contents with project root as root folder."""

    here = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(here, name)
    with io.open(path, "r", encoding="utf-8") as fd:
        content = fd.read()
    if splitlines:
        content = [row for row in content.splitlines() if row]
    return content


# Define some helper lists.
basenames = [
    "countries",
    "countriesmeta",
    "gshhs",
    "gshhsmeta",
    "rivers",
    "riversmeta",
    "states",
    "statesmeta"
]
resolutions = [
    "c",
    "l",
    "i",
    "h",
    "f",
]
grids = [
    "1.25",
    "2.5",
    "5",
    "10",
]

# Define data assets.
data_dat_files = [
    "%s_%s.dat" % (basename, res)
    for basename, res in itertools.product(basenames, resolutions[:3])
]
data_bin_files = [
    "lsmask_%smin_%s.bin" % (grid, res)
    for grid, res in itertools.product(grids, resolutions)
]
data_usc_files = [
    "UScounties.%s" % ext
    for ext in ("dbf", "prj", "shp", "shx")
]
data_other_files = [
    "epsg",
    "bmng.jpg",
    "etopo1.jpg",
    "shadedrelief.jpg",
]

data_files = data_dat_files + data_bin_files + data_usc_files + data_other_files

setup(**{
    "name":
        "basemap_data",
    "version":
        "1.3.2",
    "license":
        "GNU Lesser General Public License v3 or later (LGPLv3+)",
    "description":
        "Data assets for matplotlib basemap",
    "long_description":
        get_content("README.md"),
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
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
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
        "mpl_toolkits.basemap_data",
    ],
    "package_dir":
        {"": "src"},
    "packages":
        find_packages(where="src"),
    "package_data": {
        "mpl_toolkits.basemap_data":
            data_files,
    },
    "python_requires":
        ", ".join([
            ">=2.6",
            "!=3.0.*",
            "!=3.1.*",
            "<4",
        ]),
    "project_urls": {
        "Bug Tracker":
            "https://github.com/matplotlib/basemap/issues",
        "Documentation":
            "https://matplotlib.org/basemap/",
        "Source":
            "https://github.com/matplotlib/basemap",
    },
})
