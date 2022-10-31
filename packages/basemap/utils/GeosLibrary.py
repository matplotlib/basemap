#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Víctor Molina García

# GeosLibrary.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GeosLibrary.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with GeosLibrary.py. If not, see <https://www.gnu.org/licenses/>.
#
"""Helper class to download and build the GEOS library."""

import io
import os
import sys
import ssl
import glob
import shutil
import struct
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

        self.version_tuple = tuple(map(int, version.split(".")))

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

    @property
    def version(self):
        """GEOS library version in string format."""

        return ".".join(map(str, self.version_tuple))

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

        with contextlib.closing(conn):
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

    def extract(self, overwrite=True):
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
            shutil.rmtree(zipfold)

        # Decompress zip file.
        with contextlib.closing(ZipFile(zippath, "r")) as fd:
            fd.extractall(self.root)

        # Ensure that GEOS internal sh scripts can be executed.
        for path in sorted(glob.glob(os.path.join(zipfold, "tools", "*.sh"))):
            os.chmod(path, 0o755)

        # Apply specific patches for GEOS < 3.6.0.
        if self.version_tuple < (3, 6, 0):
            # The SVN revision file is not created on the fly before 3.6.0.
            svn_hfile = os.path.join(zipfold, "geos_svn_revision.h")
            if not os.path.exists(svn_hfile):
                with io.open(svn_hfile, "wb") as fd:
                    text = "#define GEOS_SVN_REVISION 0"
                    fd.write(text.encode())
            # Reduce warnings when compiling with `nmake` on Windows.
            cmakefile = os.path.join(zipfold, "CMakeLists.txt")
            if os.path.exists(cmakefile):
                with io.open(cmakefile, "r", encoding="utf-8") as fd:
                    lines = fd.readlines()
                with io.open(cmakefile, "wb") as fd:
                    oldtext = 'string(REGEX REPLACE "/W[0-9]" "/W4"'
                    newtext = oldtext.replace("W4", "W1")
                    for line in lines:
                        fd.write(line.replace(oldtext, newtext).encode())

        # Apply specific patches for 3.6.0 <= GEOS < 3.7.0 on Windows.
        if (3, 6, 0) <= self.version_tuple < (3, 7, 0) and os.name == "nt":
            autogen_file = os.path.join(zipfold, "autogen.bat")
            subprocess.call([autogen_file], cwd=zipfold)
            cppfile = os.path.join(zipfold, "src", "geomgraph", "DirectedEdgeStar.cpp")
            with io.open(cppfile, "r", encoding="utf-8") as fd:
                lines = fd.readlines()
            with io.open(cppfile, "wb") as fd:
                oldtext = "DirectedEdgeStar::print() const"
                newtext = oldtext.replace(" const", "")
                for line in lines:
                    fd.write(line.replace(oldtext, newtext).encode())
            hfile = os.path.join(zipfold, "include", "geos", "geomgraph", "DirectedEdgeStar.h")
            with io.open(hfile, "r", encoding="utf-8") as fd:
                lines = fd.readlines()
            with io.open(hfile, "wb") as fd:
                oldtext = "virtual std::string print() const;"
                newtext = oldtext.replace(" const", "")
                for line in lines:
                    fd.write(line.replace(oldtext, newtext).encode())

        # Patch CMakeLists to link shared geos_c with static geos.
        if self.version_tuple < (3, 8, 0):
            cmakefile = os.path.join(zipfold, "capi", "CMakeLists.txt")
            oldtext = "target_link_libraries(geos_c geos)"
            newtext = "target_link_libraries(geos_c geos-static)"
        else:
            cmakefile = os.path.join(zipfold, "CMakeLists.txt")
            oldtext = 'add_library(geos "")'
            newtext = 'add_library(geos STATIC "")'
        with io.open(cmakefile, "r", encoding="utf-8") as fd:
            lines = fd.readlines()
        with io.open(cmakefile, "wb") as fd:
            found_sharedline = False
            shared_oldtext = "if(BUILD_SHARED_LIBS)"
            shared_newtext = "if(FALSE)"
            for line in lines:
                if not found_sharedline and shared_oldtext in line:
                    line = line.replace(shared_oldtext, shared_newtext)
                    found_sharedline = True
                fd.write(line.replace(oldtext, newtext).encode())

        # Patch doc CMakeLists in GEOS 3.8.x series.
        if (3, 8, 0) <= self.version_tuple < (3, 9, 0):
            cmakefile = os.path.join(zipfold, "doc", "CMakeLists.txt")
            oldtext1 = "target_include_directories(test_geos_unit\n"
            newtext1 = "if(BUILD_TESTING)\n    {0}".format(oldtext1)
            oldtext2 = "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)\n"
            newtext2 = "{0}endif()\n".format(oldtext2)
            with io.open(cmakefile, "r", encoding="utf-8") as fd:
                lines = fd.readlines()
            with io.open(cmakefile, "wb") as fd:
                for line in lines:
                    line = line.replace(oldtext1, newtext1)
                    line = line.replace(oldtext2, newtext2)
                    fd.write(line.encode())

    def build(self, installdir=None, toolset=None, njobs=1):
        """Build and install GEOS from source."""

        # Download and extract zip file if not present.
        zipfold = os.path.join(self.root, "geos-{0}".format(self.version))
        self.extract(overwrite=True)
        version = self.version_tuple

        # Define build and install directory.
        builddir = os.path.join(zipfold, "build")
        if installdir is None:
            installdir = os.path.expanduser("~/.local/share/libgeos")
        installdir = os.path.abspath(installdir)

        # Define generic configure and build options.
        config_opts = [
            "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_INSTALL_PREFIX={0}".format(installdir),
            "-D{0}=OFF".format("GEOS_ENABLE_TESTS" if version < (3, 8, 0)
                               else "BUILD_TESTING")
        ]
        build_opts = [
            "--config", "Release",
            "--target", "install",
        ]
        build_env = os.environ.copy()

        # Define custom configure and build options.
        if os.name == "nt":
            config_opts += ["-DCMAKE_CXX_FLAGS='/wd4251 /wd4458 /wd4530 /EHsc'"]
            if version >= (3, 6, 0) and sys.version_info[:2] >= (3, 3):
                if toolset is not None:
                    config_opts += ["-DCMAKE_GENERATOR_TOOLSET={0}".format(toolset)]
                build_opts = ["-j", "{0:d}".format(njobs)] + build_opts
            else:
                win64 = (8 * struct.calcsize("P") == 64)
                config_opts = ["-G", "NMake Makefiles"] + config_opts
                build_opts.extend([
                    "--",
                    "WIN64={0}".format("YES" if win64 else "NO"),
                    "BUILD_BATCH={0}".format("YES" if njobs > 1 else "NO"),
                ])
                if sys.version_info[:2] < (3, 3):
                    build_opts += ["MSVC_VER=1500"]
        else:
            build_env["MAKEFLAGS"] = "-j {0:d}".format(njobs)
            if version >= (3, 7, 0):
                config_opts += ["-DCMAKE_CXX_FLAGS='-fPIC'"]

        # Call cmake configure after ensuring that the build directory exists.
        try:
            os.makedirs(builddir)
        except OSError:
            pass
        subprocess.call(["cmake", ".."] + config_opts, cwd=builddir)

        # Call cmake build after ensuring that the install directory exists.
        try:
            os.makedirs(installdir)
        except OSError:
            pass
        subprocess.call(["cmake", "--build", "."] + build_opts,
                        cwd=builddir, env=build_env)
