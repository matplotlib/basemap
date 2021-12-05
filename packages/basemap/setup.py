#! /usr/bin/env python
# -*- coding: utf8 -*-
# flake8: noqa: E122
"""basemap -- Plot data on map projections with matplotlib."""

import io
import os
import sys
import glob
import warnings
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


def checkversion(directory):
    """Return GEOS version from GEOS C-API header file (geos_c.h)."""

    version = None
    try:
        header_path = os.path.join(directory, "include", "geos_c.h")
        with io.open(header_path, "r", encoding="utf-8") as fd:
            for line in fd:
                if line.startswith("#define GEOS_VERSION"):
                    version = line.split()[2]
                    break
    except IOError:
        pass
    return version


# Initialise include and library dirs.
include_dirs = []
library_dirs = []

# Define NumPy include dirs.
numpy_include_path = os.environ.get("NUMPY_INCLUDE_PATH", None)
if numpy_include_path is not None:
    include_dirs.append(numpy_include_path)
else:
    try:
        import numpy
        include_dirs.append(numpy.get_include())
    except ImportError as err:
        warnings.warn("unable to locate NumPy headers", RuntimeWarning)

# Define GEOS install directory (from environment variable or trying to guess).
geos_installdir = os.environ.get("GEOS_DIR", None)
if geos_installdir is None:
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
        geos_version = checkversion(folder)
        if geos_version is not None and geos_version >= '"3.1.1"':
            geos_installdir = folder
            break

# Define GEOS include, library and runtime dirs.
if geos_installdir is None:
    warnings.warn(" ".join([
        "Cannot find GEOS library in standard locations ('{0}').",
        "Please install the corresponding packages using your",
        "software management system or set the environment variable",
        "GEOS_DIR to point to the location where GEOS is installed",
        "(for example, if 'geos_c.h' is in '/usr/local/include'",
        "and 'libgeos_c' is in '/usr/local/lib', then you need to",
        "set GEOS_DIR to '/usr/local'",
    ]).format("', '".join(geos_search_locations)), RuntimeWarning)
else:
    include_dirs.append(os.path.join(geos_installdir, "include"))
    library_dirs.append(os.path.join(geos_installdir, "lib"))
    library_dirs.append(os.path.join(geos_installdir, "lib64"))
    runtime_library_dirs = library_dirs
    data_files = []
    if os.name == "nt":
        # On Windows:
        # - DLLs get installed under `bin`.
        # - We need to inject the DLL in the wheel using `data_files`.
        # - We do not use `runtime_library_dirs` as workaround for a
        #   `distutils` bug (http://bugs.python.org/issue2437).
        library_dirs.append(os.path.join(geos_installdir, "bin"))
        runtime_library_dirs = []
        dlls = glob.glob(os.path.join(geos_installdir, "*", "geos_c.dll"))
        if dlls:
            data_files.append(("../..", sorted(dlls)))

# Define `_geoslib` extension module. It cannot be installed in the
# `mpl_toolkits.basemap` namespace or `Basemap` objects will not be pickleable.
ext_modules = [
    Extension(**{
        "name":
            "_geoslib",
        "sources": [
            "src/_geoslib.pyx",
        ],
        "libraries": [
            "geos_c",
        ],
        "include_dirs":
            include_dirs,
        "library_dirs":
            library_dirs,
        "runtime_library_dirs":
            runtime_library_dirs,
    }),
]

# Define all the different requirements.
dev_requires = get_content("requirements-dev.txt", splitlines=True)
doc_requires = get_content("requirements-doc.txt", splitlines=True)
setup_requires = get_content("requirements-setup.txt", splitlines=True)
install_requires = get_content("requirements.txt", splitlines=True)
if sys.version_info[:2] == (3, 2):
    # Hack for Python 3.2 because pip < 8 cannot handle version markers.
    marker = '; python_version == "3.2"'
    dev_requires = [
        item.replace(marker, "") for item in dev_requires
        if item.endswith(marker) or "python_version" not in item]
    doc_requires = [
        item.replace(marker, "") for item in doc_requires
        if item.endswith(marker) or "python_version" not in item]
    setup_requires = [
        item.replace(marker, "") for item in setup_requires
        if item.endswith(marker) or "python_version" not in item]
    install_requires = [
        item.replace(marker, "") for item in install_requires
        if item.endswith(marker) or "python_version" not in item]

# To create the source .tar.gz file:  python setup.py sdist
# To create the universal wheel file: python setup.py bdist_wheel --universal
setup(**{
    "name":
        "basemap",
    "version":
        "1.3.0b1",
    "license":
        "MIT",
    "description":
        "Plot data on map projections with matplotlib",
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
        "License :: OSI Approved :: MIT License",
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
        {"": "src"},
    "packages":
        find_packages(where="src"),
    "ext_modules":
        ext_modules,
    "data_files":
        data_files,
    "python_requires":
        ", ".join([
            ">=2.6",
            "!=3.0.*",
            "!=3.1.*",
            "<4",
        ]),
    "setup_requires":
        setup_requires,
    "install_requires":
        install_requires,
    "extras_require": {
        "dev":
            dev_requires,
        "doc":
            doc_requires,
    },
    "project_urls": {
        "Bug Tracker":
            "https://github.com/matplotlib/basemap/issues",
        "Documentation":
            "https://matplotlib.org/basemap/",
        "Source":
            "https://github.com/matplotlib/basemap",
    },
})
