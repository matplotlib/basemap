#! /usr/bin/env python
# -*- coding: utf-8 -*-
# flake8: noqa: E122
"""basemap -- Plot data on map projections with matplotlib."""

import io
import os
import sys
import glob
import warnings
from setuptools import setup
from setuptools import find_packages
from setuptools.command.sdist import sdist
from setuptools.dist import Distribution
from setuptools.extension import Extension


def get_content(name, splitlines=False):
    """Return the file contents with project root as root folder."""

    here = os.path.abspath(os.path.dirname(__file__))
    path = os.path.join(here, name)
    with io.open(path, "r", encoding="utf-8") as fd:
        content = fd.read()
    if splitlines:
        content = [row for row in content.splitlines() if row]
    return content


def get_geos_install_prefix():
    """Return GEOS installation prefix or None if not found."""

    env_candidate = os.environ.get("GEOS_DIR", None)
    if env_candidate is not None:
        candidates = [env_candidate]
    else:
        candidates = [os.path.expanduser("~/local"), os.path.expanduser("~"),
                      "/usr/local", "/usr", "/opt/local", "/opt", "/sw"]

    # Prepare filename pattern to find the GEOS library.
    extensions = {"win32": "dll", "cygwin": "dll", "darwin": "dylib"}
    libext = extensions.get(sys.platform, "so*")
    libname = "*geos_c*.{0}".format(libext)
    libdirs = ["bin", "lib", "lib/x86_64-linux-gnu", "lib64"]

    for prefix in candidates:
        libfiles = []
        for libdir in libdirs:
            libfiles.extend(glob.glob(os.path.join(prefix, libdir, libname)))
        hfile = os.path.join(prefix, "include", "geos_c.h")
        if os.path.isfile(hfile) and libfiles:
            return prefix

    # At this point, the GEOS library was not found, so we throw a warning if
    # the user is trying to build the library.
    build_cmds = ("bdist_wheel", "build", "install")
    if any(cmd in sys.argv[1:] for cmd in build_cmds):
        warnings.warn(" ".join([
            "Cannot find GEOS library and/or headers in standard locations",
            "('{0}'). Please install the corresponding packages using your",
            "software management system or set the environment variable",
            "GEOS_DIR to point to the location where GEOS is installed",
            "(for example, if 'geos_c.h' is in '/usr/local/include'",
            "and 'libgeos_c' is in '/usr/local/lib', then you need to",
            "set GEOS_DIR to '/usr/local'",
        ]).format("', '".join(candidates)), RuntimeWarning)
    return None


class basemap_sdist(sdist):
    """Custom `sdist` so that it will not pack DLLs on Windows if present."""

    def run(self):
        """Custom `run` command."""

        # Replace DLL data files and add GEOS build script.
        orig_data_files = self.distribution.data_files
        self.distribution.data_files = [
            (".", glob.glob(os.path.join("utils", "*.py")))]

        # Run the original `run` method and leave `data_files` as it was found.
        try:
            sdist.run(self)
        finally:
            self.distribution.data_files = orig_data_files


# Initialise include and library dirs.
data_files = []
include_dirs = []
library_dirs = []
runtime_library_dirs = []

# Define NumPy include dirs.
numpy_include_path = os.environ.get("NUMPY_INCLUDE_PATH", None)
if numpy_include_path is not None:
    include_dirs.append(numpy_include_path)
else:
    try:
        import numpy
        include_dirs.append(numpy.get_include())
    except ImportError as err:
        build_cmds = ("bdist_wheel", "build", "install")
        if any(cmd in sys.argv[1:] for cmd in build_cmds):
            warnings.warn("unable to locate NumPy headers", RuntimeWarning)

# Define GEOS include, library and runtime dirs.
geos_install_prefix = get_geos_install_prefix()
if geos_install_prefix is not None:
    include_dirs.append(os.path.join(geos_install_prefix, "include"))
    library_dirs.append(os.path.join(geos_install_prefix, "lib"))
    library_dirs.append(os.path.join(geos_install_prefix, "lib64"))
    runtime_library_dirs = library_dirs
    if os.name == "nt":
        # On Windows:
        # - DLLs get installed under `bin`.
        # - We need to inject later the DLL in the wheel using `data_files`.
        # - We do not use `runtime_library_dirs` as workaround for a
        #   `distutils` bug (http://bugs.python.org/issue2437).
        library_dirs.append(os.path.join(geos_install_prefix, "bin"))
        runtime_library_dirs = []
        dlls = glob.glob(os.path.join(geos_install_prefix, "*", "geos_c.dll"))
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
for ext in ext_modules:
    ext.cython_directives = [
        ("language_level", str(sys.version_info[0])),
    ]

# Define all the different requirements.
setup_requires = get_content("requirements-setup.txt", splitlines=True)
install_requires = get_content("requirements.txt", splitlines=True)
if sys.version_info[:2] == (3, 2):
    # Hack for Python 3.2 because pip < 8 cannot handle version markers.
    marker1 = '; python_version == "3.2"'
    marker2 = '; python_version >= "2.7"'
    setup_requires = [
        item.replace(marker1, "").replace(marker2, "") for item in setup_requires
        if item.endswith(marker1) or item.endswith(marker2)
        or "python_version" not in item]
    install_requires = [
        item.replace(marker1, "").replace(marker2, "") for item in install_requires
        if item.endswith(marker1) or item.endswith(marker2)
        or "python_version" not in item]
else:
    marker1 = '; python_version == "3.2"'
    setup_requires = [item for item in setup_requires if not item.endswith(marker1)]
    install_requires = [item for item in install_requires if not item.endswith(marker1)]

setup(**{
    "name":
        "basemap",
    "version":
        "1.3.7",
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
        "doc":
            get_content("requirements-doc.txt", splitlines=True),
        "lint":
            get_content("requirements-lint.txt", splitlines=True),
        "test":
            get_content("requirements-test.txt", splitlines=True),
    },
    "cmdclass": {
        "sdist": basemap_sdist,
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
