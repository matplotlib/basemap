#! /usr/bin/env python
"""Helper script to be run by `cibuildwheel` as `before_all` step."""

import os
import sys

HERE = os.path.abspath(__file__)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "packages", "basemap"))
import utils  # noqa: E402  # pylint: disable=imports


def main():
    """Build the GEOS library based on parsed environment variables."""

    geos_version = os.environ.get("GEOS_VERSION", None)
    if geos_version is None:
        raise ValueError("Undefined environment variable GEOS_VERSION")

    geos_dir = os.environ.get("GEOS_DIR", None)
    if geos_dir is None:
        raise ValueError("Undefined environment variable GEOS_DIR")

    geos_njobs = int(os.environ.get("GEOS_NJOBS", 1))

    # pylint: disable=consider-using-f-string
    print("Running before_all script with the following settings:")
    print("GEOS_DIR: {0}".format(geos_dir))
    print("GEOS_VERSION: {0}".format(geos_version))
    print("GEOS_NJOBS: {0}".format(geos_njobs))

    utils.GeosLibrary(geos_version).build(geos_dir, njobs=geos_njobs)
    return 0


if __name__ == "__main__":
    sys.exit(main())
