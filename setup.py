import sys, glob, os, subprocess

major, minor1, minor2, s, tmp = sys.version_info
if major==2 and minor1<4 or major<2:
    raise SystemExit("""matplotlib and the basemap toolkit require Python 2.4 or later.""")

from distutils.dist import Distribution
from distutils.util import convert_path
from distutils import ccompiler, sysconfig

# Do not require numpy for just querying the package
# Taken from the netcdf-python setup file (which took it from h5py setup file).
inc_dirs = []
if any('--' + opt in sys.argv for opt in Distribution.display_option_names +
       ['help-commands', 'help']) or sys.argv[1] == 'egg_info':
    from distutils.core import setup, Extension
else:
    import numpy
    # Use numpy versions if they are available.
    from numpy.distutils.core import setup, Extension
    # append numpy include dir.
    inc_dirs.append(numpy.get_include())

try:
    from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
    from distutils.command.build_py import build_py

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
if 'GEOS_DIR' in os.environ:
    GEOS_dir = os.environ.get('GEOS_DIR')
else:
# set GEOS_dir manually here if automatic detection fails.
    GEOS_dir = None

print '1: GEOS_dir = ', GEOS_dir

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
            sys.stdout.write('GEOS lib (version %s) found in %s\n' %\
                    (geos_version[1:-1],direc))
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
    geos_include_dirs=[os.path.join(GEOS_dir,'include'),inc_dirs]
    geos_library_dirs=[os.path.join(GEOS_dir,'lib'),os.path.join(GEOS_dir,'lib64')]

# get location of pyshp lib from environment variable if it is set.
if 'PYSHP_DIR' in os.environ:
    PYSHP_dir = os.environ.get('PYSHP_DIR')
else:
    # set PYSHP_DIR to None if not defined by the user
    PYSHP_dir = None

if PYSHP_dir:
    if not PYSHP_dir in sys.path:
        sys.path.append(PYSHP_dir)
    try:
        import shapefile
    except:
        raise SystemExit("""
Could not import shapefile!'
Please make sure the pyshp module is installed and set the PYSHP_DIR
environment variable to point to the proper location.
Alternatively, you may choose to not define the PYSHP_DIR environment
variable, in which case the bundled copy of pyshp will be used.
""")
else:
    pass

print '2: GEOS_dir = ', GEOS_dir
print 'PYSHP_DIR = ', PYSHP_dir
#sys.exit(1)

# proj4 and geos extensions.
deps = glob.glob('src/*.c')
deps.remove(os.path.join('src','_proj.c'))
deps.remove(os.path.join('src','_geoslib.c'))

packages          = ['mpl_toolkits','mpl_toolkits.basemap']
namespace_packages = ['mpl_toolkits']
package_dirs       = {'':'lib'}

if PYSHP_dir is None:
    # copy shapefile into basemap dir
    import shutil
    shutil.copy('opt_lib/mpl_toolkits/basemap/shapefile.py',
                'lib/mpl_toolkits/basemap/shapefile.py')
    
extensions = [Extension("mpl_toolkits.basemap._proj",deps+['src/_proj.c'],include_dirs = ['src'],)]
# can't install _geoslib in mpl_toolkits.basemap namespace,
# or Basemap objects won't be pickleable.
if sys.platform == 'win32':
# don't use runtime_library_dirs on windows (workaround
# for a distutils bug - http://bugs.python.org/issue2437).
    #extensions.append(Extension("mpl_toolkits.basemap._geoslib",['src/_geoslib.c'],
    extensions.append(Extension("_geoslib",['src/_geoslib.c'],
                                library_dirs=geos_library_dirs,
                                include_dirs=geos_include_dirs,
                                libraries=['geos']))
else:
    #extensions.append(Extension("mpl_toolkits.basemap._geoslib",['src/_geoslib.c'],
    extensions.append(Extension("_geoslib",['src/_geoslib.c'],
                                library_dirs=geos_library_dirs,
                                runtime_library_dirs=geos_library_dirs,
                                include_dirs=geos_include_dirs,
                                libraries=['geos_c']))

# Specify all the required mpl data
# create pyproj binary datum shift grid files.
pathout =\
os.path.join('lib',os.path.join('mpl_toolkits',os.path.join('basemap','data')))
if sys.argv[1] not in ['sdist','clean']:
    cc = ccompiler.new_compiler()
    sysconfig.get_config_vars()
    sysconfig.customize_compiler(cc)
    cc.set_include_dirs(['src'])
    objects = cc.compile(['nad2bin.c', 'src/pj_malloc.c'])
    execname = 'nad2bin'
    cc.link_executable(objects, execname)
    llafiles = glob.glob('datumgrid/*.lla')
    cmd = os.path.join(os.getcwd(),execname)
    for f in llafiles:
        fout = os.path.basename(f.split('.lla')[0])
        fout = os.path.join(pathout,fout)
        strg = '%s %s < %s' % (cmd, fout, f)
        sys.stdout.write('executing %s\n' % strg)
        subprocess.call(strg,shell=True)
datafiles = glob.glob(os.path.join(pathout,'*'))
datafiles = [os.path.join('data',os.path.basename(f)) for f in datafiles]
package_data = {'mpl_toolkits.basemap':datafiles}

__version__ = "1.0.8"
setup(
  name              = "basemap",
  version           = __version__,
  description       = "Plot data on map projections with matplotlib",
  long_description  = """
  An add-on toolkit for matplotlib that lets you plot data
  on map projections with coastlines, lakes, rivers and political boundaries.
  See http://www.scipy.org/wikis/topical_software/Maps for an
  example of what it can do.""",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  download_url      = "https://downloads.sourceforge.net/project/matplotlib/matplotlib-toolkits/basemap-{0}/basemap-{0}.tar.gz".format(__version__),
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  install_requires  = ["numpy>=1.2.1", "matplotlib>=1.0.0"],
  platforms         = ["any"],
  license           = "OSI Approved",
  keywords          = ["python","plotting","plots","graphs","charts","GIS","mapping","map projections","maps"],
  classifiers       = ["Development Status :: 5 - Production/Stable",
                       "Intended Audience :: Science/Research",
                       "License :: OSI Approved",
                       "Programming Language :: Python",
                       "Programming Language :: Python :: 3",
                       "Topic :: Scientific/Engineering :: Visualization",
                       "Topic :: Software Development :: Libraries :: Python Modules",
                       "Operating System :: OS Independent"],
  packages          = packages,
  namespace_packages = namespace_packages,
  package_dir       = package_dirs,
  ext_modules       = extensions,
  cmdclass = {'build_py': build_py},
  package_data = package_data
  )

# clean up afterwards
if PYSHP_dir is None:
    os.remove('lib/mpl_toolkits/basemap/shapefile.py')
