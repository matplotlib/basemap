#! /usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa: E122
"""basemap_data_hires -- High-resolution data assets for matplotlib basemap."""

import io
import os
import itertools
from setuptools import setup
from setuptools import find_namespace_packages


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

# Define data assets.
data_files = [
    "%s_%s.dat" % (basename, res)
    for basename, res in itertools.product(basenames, resolutions[3:])
]

setup(**{
    "name":
        "basemap_data_hires",
    "version":
        "2.0.0-dev",
    "license":
        "GNU Lesser General Public License v3 or later (LGPLv3+)",
    "description":
        "High-resolution data assets for matplotlib basemap",
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
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    "keywords": [
        "GIS",
        "maps",
        "plots",
    ],
    "package_dir":
        {"": "src"},
    "packages":
        find_namespace_packages(where="src"),
    "package_data": {
        "mpl_toolkits.basemap_data":
            data_files,
    },
    "python_requires":
        ", ".join([
            ">=3.9",
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
