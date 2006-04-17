from distutils.core import setup, Extension
from distutils.util import convert_path
import sys, glob, os

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

extensions = [Extension("proj4",deps,include_dirs = ['src'],)]
packages          = ['matplotlib.toolkits.basemap']
package_dirs       = {'':'lib'}

# only install shapelib and dbflib if user hasn't got them.
#try: import shapelib
#except ImportError: haveshapelib = False
#else: haveshapelib = True
#try: import dbflib
#except ImportError: havedbflib = False
#else: havedbflib = True
#if not haveshapelib or not havedbflib:
# always intall shapelib and dbflib.
if 1:
    packages = packages + ['shapelib','dbflib']
    package_dirs['shapelib'] ='lib/shapelib'
    package_dirs['dbflib'] ='lib/dbflib'
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
else:
    additional_params = {}

datadir = os.environ.get('BASEMAP_DATA_PATH')
if not datadir:
    datadir ='share/basemap'

setup(
  name              = "basemap",
  version           = "0.9",
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
  data_files        = [(datadir,['data/countries_c.txt','data/states_c.txt','data/rivers_c.txt','data/gshhs_c.txt','data/countries_l.txt','data/states_l.txt','data/rivers_l.txt','data/gshhs_l.txt'])],
  **additional_params
  )
