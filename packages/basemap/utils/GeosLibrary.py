#! /usr/bin/env python
"""Helper class to download and build the GEOS library."""

import io
import os
import ssl
import glob
import shutil
import tempfile
import contextlib
import subprocess
import datetime as dt
from zipfile import ZipFile
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen


URL_DATETIME_FMT = "%a, %d %b %Y %H:%M:%S GMT"
GEOS_BASEURL = "https://github.com/libgeos/geos/archive/refs/tags"


class GeosLibrary(object):
    """Helper class to download, build and install GEOS."""

    def __init__(self, version, root=None):
        """Initialise a new :class:`GeosLibrary` instance."""

        self.version = version
        if root is None:
            self.temp = True
            self.root = tempfile.mkdtemp(prefix="tmp_geoslibrary_")
        else:
            self.temp = False
            self.root = root
            try:
                os.makedirs(self.root)
            except OSError:
                pass

    def __del__(self):
        """Clean up after :class:`GeosLibrary` destruction."""

        if self.temp:
            try:
                shutil.rmtree(self.root)
            except OSError:
                pass

    def download(self):
        """Download GEOS zip source code into :class:`GeosLibrary` root."""

        # Define download link.
        link = "{0}/{1}.zip".format(GEOS_BASEURL, self.version)
        suffix = os.path.splitext(link)[-1]

        # Define output path.
        zipname = "geos-{0}".format(link.rsplit("/", 1)[-1])
        zippath = os.path.join(self.root, zipname)

        # Handle creation of the HTTP request.
        kwargs = {}
        if hasattr(ssl, "SSLContext") and hasattr(ssl, "PROTOCOL_TLSv1_2"):
            kwargs.update(context=ssl.SSLContext(ssl.PROTOCOL_TLSv1_2))
        try:
            conn = urlopen(link, **kwargs)
        except TypeError:
            # Fallback if `urlopen` does not accept context.
            conn = urlopen(link)

        with conn:
            # Try to get the file timestamp from the HTTP request header.
            date = conn.headers.get("Last-Modified")
            if date is not None:
                date = dt.datetime.strptime(date, URL_DATETIME_FMT)
                date = date.replace(tzinfo=dt.timezone.utc)
                date = date.timestamp()
            with tempfile.NamedTemporaryFile(suffix=suffix) as tmpzipobj:
                # Copy the buffer into a temporary file.
                shutil.copyfileobj(conn, tmpzipobj)
                # Move the file descriptor pointer to the beginning and
                # simply copy it to the final destination.
                tmpzipobj.seek(0)
                with open(zippath, "wb") as zipobj:
                    shutil.copyfileobj(tmpzipobj, zipobj)
            # Assign the timestamps to the final file if available.
            if date is not None:
                os.utime(zippath, (date, date))

    def decompress(self, overwrite=True):
        """Decompress GEOS zip source code into :class:`GeosLibrary` root."""

        # Download zip file if not present.
        zippath = os.path.join(self.root, "geos-{0}.zip".format(self.version))
        if not os.path.exists(zippath):
            self.download()

        # Remove destination folder if present and requested.
        zipfold = os.path.join(self.root, "geos-{0}".format(self.version))
        if os.path.exists(zipfold):
            if not overwrite:
                raise OSError("folder '{0}' already exists".format(zipfold))
            else:
                shutil.rmtree(zipfold)

        # Decompress zip file.
        with contextlib.closing(ZipFile(zippath, "r")) as fd:
            fd.extractall(self.root)

        # Ensure that GEOS internal sh scripts can be executed.
        for path in sorted(glob.glob(os.path.join(zipfold, "tools", "*.sh"))):
            os.chmod(path, 0o755)

        # Patch CMakeLists so that libgeos_c.so does not depend on libgeos.so.
        cmakefile = os.path.join(zipfold, "capi", "CMakeLists.txt")
        with io.open(cmakefile, "r", encoding="utf-8") as fd:
            lines = fd.readlines()
        with io.open(cmakefile, "wb") as fd:
            oldtext = "target_link_libraries(geos_c geos)"
            newtext = "target_link_libraries(geos_c geos-static)"
            for line in lines:
                fd.write(line.replace(oldtext, newtext).encode())

    def build(self, njobs=1, installdir=None):
        """Build and install GEOS from source."""

        # Download and decompress zip file if not present.
        zipfold = os.path.join(self.root, "geos-{0}".format(self.version))
        self.decompress(overwrite=True)

        # Define installation directory.
        if installdir is None:
            installdir = os.path.expanduser("~/.local/share/libgeos")
        installdir = os.path.abspath(installdir)

        # Create installation directory if not present.
        try:
            os.makedirs(installdir)
        except OSError:
            pass

        # Now move to the GEOS source code folder and build with CMake.
        cwd = os.getcwd()
        try:
            os.chdir(zipfold)
            subprocess.call(["cmake",
                             "-S", ".",
                             "-B", "build",
                             "-DCMAKE_BUILD_TYPE=Release",
                             "-DCMAKE_INSTALL_PREFIX={0}".format(installdir)])
            subprocess.call(["cmake",
                             "--build", "build", "-j", str(int(njobs)),
                             "--target", "install"])
        finally:
            os.chdir(cwd)
