# basemap build options can be modified with the setup.cfg file. See
# setup.cfg for more information.
import sys, glob, os, numpy
major, minor1, minor2, s, tmp = sys.version_info
if major==2 and minor1<=3:
    # setuptools monkeypatches distutils.core.Distribution to support
    # package_data
    try: import setuptools
    except ImportError:
        raise SystemExit("""
matplotlib requires setuptools for installation.  Please download
http://peak.telecommunity.com/dist/ez_setup.py and run it (as su if
you are doing a system wide install) to install the proper version of
setuptools for your system""")
#from distutils.core import setup
#from distutils.core import Extension
from numpy.distutils.core  import setup, Extension
from distutils.util import convert_path

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

def checkversion(GEOS_dir):
    """check geos C-API header file (geos_c.h)"""
    try:
        f = open(os.path.join(GEOS_dir,'include/geos_c.h'))
    except IOError:
        return None
    geos_version = None
    for line in f:
        if line.startswith('#define GEOS_VERSION'):
            geos_version = line.split()[2]
    return geos_version

# get location of geos lib from environment variable if it is set.
if os.environ.has_key('GEOS_DIR'):
    GEOS_dir = os.environ.get('GEOS_DIR')
else:
# set GEOS_dir manually here if automatic detection fails.
    GEOS_dir = None

if GEOS_dir is None:
    # if GEOS_dir not set, check a few standard locations.
    GEOS_dirs = ['/usr/local','/sw','/opt','/opt/local',os.path.expanduser('~')]
    for direc in GEOS_dirs:
        geos_version = checkversion(direc)
        print 'checking for GEOS lib in %s ....' % direc
        if geos_version != '"2.2.3"':
            continue
        else:
            print 'GEOS lib found in %s' % direc
            GEOS_dir = direc
            break
else:
    geos_version = checkversion(GEOS_dir)

if geos_version != '"2.2.3"' or GEOS_dir is None:
    raise SystemExit("""
Can't find geos library version 2.2.3. Please set the
environment variable GEOS_DIR to point to the location
where geos 2.2.3 is installed (for example, if geos_c.h
is in /usr/local/include, and libgeos_c is in /usr/local/lib,
set GEOS_DIR to /usr/local), or edit the setup.py script
manually and set the variable GEOS_dir (right after the line
that says "set GEOS_dir manually here".""")
else:
    geos_include_dirs=[os.path.join(GEOS_dir,'include'),numpy.get_include()]
    geos_library_dirs=[os.path.join(GEOS_dir,'lib'),os.path.join(GEOS_dir,'lib64')]

# proj4 and geos extensions.
deps = glob.glob('src/*.c')
deps.remove(os.path.join('src','_proj.c'))
deps.remove(os.path.join('src','_geod.c'))
deps.remove(os.path.join('src','_geoslib.c'))

packages          = ['mpl_toolkits','mpl_toolkits.basemap']
package_dirs       = {'':'lib'}
extensions = [Extension("mpl_toolkits.basemap._proj",deps+['src/_proj.c'],include_dirs = ['src'],)]
extensions.append(Extension("mpl_toolkits.basemap._geod",deps+['src/_geod.c'],include_dirs = ['src'],))
# for some reason, pickling won't work if this extension is installed
# as "matplotlib.toolkits.basemap._geoslib"
extensions.append(Extension("_geoslib",['src/_geoslib.c'],
                            library_dirs=geos_library_dirs,
                            runtime_library_dirs=geos_library_dirs,
                            include_dirs=geos_include_dirs,
                            libraries=['geos_c','geos']))

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

# check setup.cfg file to see how to install auxilliary packages.
options = {}
if os.path.exists("setup.cfg"):
    import ConfigParser
    config = ConfigParser.SafeConfigParser()
    config.read("setup.cfg")
    try: options['provide_pydap'] = config.getboolean("provide_packages", "pydap")
    except: options['provide_pydap'] = 'auto'
    try: options['provide_httplib2'] = config.getboolean("provide_packages", "httplib2")
    except: options['provide_httplib2'] = 'auto'
else:
    options['provide_pydap'] = 'auto'
    options['provide_httplib2'] = 'auto'

provide_pydap = options['provide_pydap']
if provide_pydap == 'auto': # install pydap stuff if not already available.
   # only the client is installed (not the server).
    __dapversion__ = None
    print 'checking to see if required version of pydap installed ..'
    try:
        from dap.lib import __version__ as __dapversion__
    except ImportError:
        print 'pydap not installed, client will be installed'
        packages = packages + ['dap','dap.util','dap.parsers']
        package_dirs['dap'] = os.path.join('lib','dap')
    else:
        print 'pydap installed, checking version ...'
    # install dap client anyway if installed version is older than
    # version provided here.
    if __dapversion__ is not None:
        __dapversion__ = [repr(v)+'.' for v in __dapversion__]
        __dapversion__ = ''.join(__dapversion__)[:-1]
        if __dapversion__ < '2.2.6.2':
            print 'required version of pydap not installed, client will be installed'
            packages = packages + ['dap','dap.util','dap.parsers']
            package_dirs['dap'] = os.path.join('lib','dap')
        else:
            print 'pydap version OK, will not be installed'
elif provide_pydap: # force install of pydap stuff.
    print 'forcing install of included pydap client'
    packages = packages + ['dap','dap.util','dap.parsers']
    package_dirs['dap'] = os.path.join('lib','dap')
else:
    print 'will not install pydap'

provide_httplib2 = options['provide_httplib2']
if provide_httplib2  == 'auto':
    print 'checking to see if httplib2 installed ..'
    try:
        import httplib2
    except ImportError:
        print 'httplib2 not installed, will be installed'
        packages = packages + ['httplib2']
        package_dirs['httlib2'] = os.path.join('lib','httplib2')
    else:
        print 'httplib2 installed'
elif provide_httplib2: # force install of httplib2
    print 'forcing install of included httplib2'
    packages = packages + ['httplib2']
    package_dirs['httlib2'] = os.path.join('lib','httplib2')
else:
    print 'will not install httplib2'

# Specify all the required mpl data
pyproj_datafiles = ['data/epsg', 'data/esri', 'data/esri.extra', 'data/GL27', 'data/nad.lst', 'data/nad27', 'data/nad83', 'data/ntv2_out.dist', 'data/other.extra', 'data/pj_out27.dist', 'data/pj_out83.dist', 'data/proj_def.dat', 'data/README', 'data/td_out.dist', 'data/test27', 'data/test83', 'data/testntv2', 'data/testvarious', 'data/world','data/bmng.jpg']
boundaryfiles = []
for resolution in ['c','l','i','h','f']:
    boundaryfiles = boundaryfiles + glob.glob("lib/mpl_toolkits/basemap/data/*_"+resolution+".dat")
boundaryfiles = [os.path.join('data',os.path.basename(bfile)) for bfile in boundaryfiles]
basemap_datafiles = boundaryfiles + ['data/5minmask.bin']
package_data = {'mpl_toolkits.basemap':pyproj_datafiles+basemap_datafiles}
setup(
  name              = "basemap",
  version           = "0.99",
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
  package_data = package_data
  )
