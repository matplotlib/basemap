import sys, glob, os
major, minor1, minor2, s, tmp = sys.version_info
if major==2 and minor1<=3:
    # setuptools monkeypatches distutils.core.Distribution to support
    # package_data
    try: import setuptools
    except ImportError:
        raise SystemExit("""\
matplotlib requires setuptools for installation.  Please download
http://peak.telecommunity.com/dist/ez_setup.py and run it (as su if
you are doing a system wide install) to install the proper version of
setuptools for your system""")
from distutils.core import Extension
from distutils.util import convert_path

import numpy


def dbf_macros():
    """Return the macros to define when compiling the dbflib wrapper.

    The returned list specifies one macro, HAVE_UPDATE_HEADER, which is
    '1' if the dbflib version we will be compiling with has the
    DBFUpdateHeader function and '0' otherwise.  To check whether
    DBFUpdateHeader is available, we scan shapefil.h for the string
    'DBFUpdateHeader'.
    """
    f = open(convert_path("pyshapelib/shapelib/shapefil.h"))
    contents = f.read()
    f.close()
    if contents.find("DBFUpdateHeader") >= 0:
        return [("HAVE_UPDATE_HEADER", "1")]
    else:
        return [("HAVE_UPDATE_HEADER", "0")]

deps = glob.glob('src/*.c')
deps.remove(os.path.join('src','_proj.c'))
deps.remove(os.path.join('src','_geod.c'))

packages          = ['matplotlib.toolkits.basemap']
package_dirs       = {'':'lib'}
extensions = [Extension("matplotlib.toolkits.basemap._proj",
                deps+['src/_proj.c'],
                include_dirs = ['src', numpy.get_include()],)]
extensions.append(Extension("matplotlib.toolkits.basemap._geod",deps+['src/_geod.c'],include_dirs = ['src'],))

# install shapelib and dbflib.
packages = packages + ['shapelib','dbflib']
package_dirs['shapelib'] = os.path.join('lib','shapelib')
package_dirs['dbflib'] = os.path.join('lib','dbflib')
extensions = extensions + \
         [Extension("shapelibc",
                    ["pyshapelib/shapelib_wrap.c",
                     "pyshapelib/shapelib/shpopen.c",
                     "pyshapelib/shapelib/shptree.c"],
                    include_dirs = ["pyshapelib/shapelib"]),
          Extension("shptree",
                    ["pyshapelib/shptreemodule.c"],
                    include_dirs = ["pyshapelib/shapelib"]),
          Extension("dbflibc",
                    ["pyshapelib/dbflib_wrap.c",
                     "pyshapelib/shapelib/dbfopen.c"],
                    include_dirs = ["pyshapelib/shapelib"],
                    define_macros = dbf_macros()) ]

if 'setuptools' in sys.modules:
# Are we running with setuptools?
# if so, need to specify all the packages in heirarchy
    additional_params = {'namespace_packages' : ['matplotlib.toolkits']}
    packages.extend(['matplotlib', 'matplotlib.toolkits'])
    setup = setuptools.setup
else:
    additional_params = {}
    from distutils.core import setup


# Specify all the required mpl data
pyproj_datafiles = ['data/epsg', 'data/esri', 'data/esri.extra', 'data/GL27', 'data/nad.lst', 'data/nad27', 'data/nad83', 'data/ntv2_out.dist', 'data/other.extra', 'data/pj_out27.dist', 'data/pj_out83.dist', 'data/proj_def.dat', 'data/README', 'data/td_out.dist', 'data/test27', 'data/test83', 'data/testntv2', 'data/testvarious', 'data/world']
basemap_datafiles = [         'data/countries_c.txt',
                              'data/states_c.txt',
                              'data/rivers_c.txt',
                              'data/gshhs_c.txt',
                              'data/countries_l.txt',
                              'data/states_l.txt',
                              'data/rivers_l.txt',
                              'data/gshhs_l.txt',
                              'data/countries_i.txt',
                              'data/states_i.txt',
                              'data/rivers_i.txt',
                              'data/gshhs_i.txt',
                              'data/5minmask.bin',
                              ]
package_data = {'matplotlib.toolkits.basemap':pyproj_datafiles+basemap_datafiles}
setup(
  name              = "basemap",
  version           = "0.9.6.1",
  description       = "Plot data on map projections with matplotlib",
  long_description  = """
  An add-on toolkit for matplotlib that lets you plot data
  on map projections with coastlines, lakes, rivers and political boundaries.
  See http://www.scipy.org/wikis/topical_software/Maps for an
  example of what it can do.""",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  download_url      = "http://sourceforge.net/projects/matplotlib",
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  platforms         = ["any"],
  license           = "OSI Approved",
  keywords          = ["python","plotting","plots","graphs","charts","GIS","mapping","map projections","maps"],
  classifiers       = ["Development Status :: 4 - Beta",
			           "Intended Audience :: Science/Research",
			           "License :: OSI Approved",
			           "Topic :: Scientific/Engineering :: Visualization",
			           "Topic :: Software Development :: Libraries :: Python Modules",
			           "Operating System :: OS Independent"],
  packages          = packages,
  package_dir       = package_dirs,
  ext_modules       = extensions,
  package_data = package_data,
  **additional_params
  )
