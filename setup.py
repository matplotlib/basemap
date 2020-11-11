#! /usr/bin/env python
# -*- coding: utf8 -*-
# flake8: noqa: E122
from __future__ import (absolute_import, division, print_function)

import re
import glob
import io
import os
import sys
from setuptools import setup
from setuptools.dist import Distribution
from setuptools.extension import Extension
__version__ = "1.3.0"

if sys.version_info < (2, 6):
    raise SystemExit("""matplotlib and the basemap toolkit require Python 2.6 or later.""")

# Do not require numpy for just querying the package
# Taken from the netcdf-python setup file (which took it from h5py setup file).
inc_dirs = []
if not (any('--' + opt in sys.argv for opt in Distribution.display_option_names +
            ['help-commands', 'help']) or sys.argv[1] == 'egg_info'):
    import numpy
    # append numpy include dir.
    inc_dirs.append(numpy.get_include())


def get_install_requirements(path):
    path = os.path.join(os.path.dirname(__file__), path)
    with io.open(path, encoding='utf-8') as fp:
        content = fp.read()
    return [req for req in content.split("\n")
            if req != '' and not req.startswith('#')]


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
# set GEOS_dir manually here if automatic detection fails.
else:
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
            sys.stdout.write('GEOS lib (version %s) found in %s\n' %
                             (geos_version[1:-1], direc))
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
    geos_include_dirs = [os.path.join(GEOS_dir, 'include')] + inc_dirs
    geos_library_dirs = [os.path.join(GEOS_dir, 'lib'), os.path.join(GEOS_dir, 'lib64')]

# can't install _geoslib in mpl_toolkits.basemap namespace,
# or Basemap objects won't be pickleable.

# don't use runtime_library_dirs on windows (workaround
# for a distutils bug - http://bugs.python.org/issue2437).
if sys.platform == 'win32':
    runtime_lib_dirs = []
else:
    runtime_lib_dirs = geos_library_dirs

extensions = [
    Extension("_geoslib", ["src/_geoslib.pyx"],
              library_dirs=geos_library_dirs,
              runtime_library_dirs=runtime_lib_dirs,
              include_dirs=geos_include_dirs,
              libraries=["geos_c"]),
]

# Define the build mode (normal, lite, data or extras).
mode_arg = [item for item in sys.argv[1:] if item.startswith("--mode")]
if len(mode_arg):
    sys.argv.remove(mode_arg[0])
    if len(mode_arg) > 1:
        raise ValueError("setup option --mode given multiple times")
mode = (mode_arg or [""])[0].lstrip("--mode").strip("=")
if mode not in ("", "lite", "data", "extras"):
    raise ValueError("invalid setup mode: {0}".format(mode))
if mode and "sdist" in sys.argv[1:]:
    raise ValueError("setup option --mode incompatible with sdist")

# Get the basemap data files.
data_pattern = os.path.join("lib", "mpl_toolkits", "basemap", "data", "*")
data_files = sorted(map(os.path.basename, glob.glob(data_pattern)))

# Define default setup parameters.
namespace_packages = [
    "mpl_toolkits",
]
packages = [
    "mpl_toolkits.basemap",
    "mpl_toolkits.basemap.data",
]
package_dirs = {"": "lib"}
package_data = {
    "mpl_toolkits.basemap.data":
        data_files,
}
install_requires = get_install_requirements("requirements.txt")

# Filter the data files depending on the mode (normal, lite, data, extras).
if mode:
    version, vbuild = (__version__.split("+") + [""])[:2]
    data_vbuild = "{0}{1}".format("+" if vbuild else "", vbuild)
    data_version = "{0}.0{1}".format(version.rsplit(".", 1)[0], data_vbuild)
    regex = re.compile("(UScounties|_[ihf]\\.dat$)")
    if mode == "lite":
        packages.remove("mpl_toolkits.basemap.data")
        package_data.pop("mpl_toolkits.basemap.data")
        install_requires.append("basemap-data == {0}".format(data_version))
    else:
        __version__ = data_version
        extensions = []
        packages.remove("mpl_toolkits.basemap")
        if mode == "data":
            data_files = [f for f in data_files if not regex.search(f)]
        elif mode == "extras":
            data_files = [f for f in data_files if regex.search(f)]
        package_data["mpl_toolkits.basemap.data"] = data_files


setup(**{
    "name":
        "basemap{0}{1}".format("-" if mode else "", mode),
    "version":
        __version__,
    "description":
        "Plot data on map projections with matplotlib",
    "long_description": """
An add-on toolkit for matplotlib that lets you plot data on map
projections with coastlines, lakes, rivers and political boundaries.
See http://matplotlib.org/basemap/users/examples.html for
examples of what it can do.""",
    "url":
        "https://matplotlib.org/basemap/",
    "download_url":
        "https://github.com/matplotlib/basemap/archive/v{0}.tar.gz".format(__version__),
    "author":
        "Jeff Whitaker",
    "author_email":
        "jeffrey.s.whitaker@noaa.gov",
    "maintainer":
        "Ben Root",
    "maintainer_email":
        "ben.v.root@gmail.com",
    "install_requires":
        install_requires,
    "platforms": [
        "any",
    ],
    "license":
        "OSI Approved",
    "keywords": [
        "python",
        "plotting",
        "plots",
        "graphs",
        "charts",
        "GIS",
        "mapping",
        "map projections",
        "maps",
    ],
    "classifiers": [
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.2",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Visualization",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    "namespace_packages":
        namespace_packages,
    "packages":
        packages,
    "package_dir":
        package_dirs,
    "package_data":
        package_data,
    "ext_modules":
        extensions,
    "python_requires":
        ", ".join([
            ">= 2.6",
            "!= 3.0.*",
            "!= 3.1.*",
            "< 4",
        ]),
})
