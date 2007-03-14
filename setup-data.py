from distutils.core import setup, Extension
import sys, os

packages          = ['matplotlib.toolkits.basemap']
package_dirs       = {'':'lib'}

if 'setuptools' in sys.modules:
    # Are we running with setuptools?
    # if so, need to specify all the packages in heirarchy
    additional_params = {'namespace_packages' : ['matplotlib.toolkits']}    
    packages.extend(['matplotlib', 'matplotlib.toolkits'])
else:
    additional_params = {}

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
    
# Specify all the required mpl data
package_data = {'matplotlib.toolkits.basemap':[
                              'data/countries_i.txt',
                              'data/states_i.txt',
                              'data/rivers_i.txt',
                              'data/gshhs_i.txt',
                              'data/countries_h.txt', 
                              'data/states_h.txt',
                              'data/rivers_h.txt',
                              'data/gshhs_h.txt',
                              ]}
setup(
  name              = "basemap-data",
  version           = "0.9.5",
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
  package_data = package_data,
  **additional_params
  )
