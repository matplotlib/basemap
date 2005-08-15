from distutils.core import setup, Extension
from distutils.util import convert_path
import sys, glob

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

extensions = [Extension("proj4",deps,include_dirs = ['src'],),
              Extension("shapelibc",
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
                        define_macros = dbf_macros())]

datadir ='share/basemap-py'+repr(sys.version_info[0])+repr(sys.version_info[1])

setup(
  name              = "basemap",
  version           = "0.6.1",
  description       = "Plot data on map projections with matplotlib",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  data_files        = [(datadir,['data/countries_c.txt','data/states_c.txt','data/countries_l.txt','data/states_l.txt','data/gshhs_c.txt','data/gshhs_l.txt','data/countries_i.txt','data/states_i.txt','data/gshhs_i.txt'])],
  packages          = ['matplotlib/toolkits','matplotlib/toolkits/basemap','shapelib','dbflib'],
  package_dir       = {'':'lib','shapelib':'pyshapelib/lib/shapelib','dbflib':'pyshapelib/lib/dbflib'},
  ext_modules       = extensions)
