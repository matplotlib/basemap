import sys, glob, os
major, minor1, minor2, s, tmp = sys.version_info
if major==2 and minor1<=3:
    # setuptools monkeypatches distutils.core.Distribution to support
    # package_data
    #try: import setuptools
    #except ImportError:
    #    raise SystemExit("""
#matplotlib requires setuptools for installation.  Please download
#http://peak.telecommunity.com/dist/ez_setup.py and run it (as su if
#you are doing a system wide install) to install the proper version of
#setuptools for your system""")
    raise SystemExit("""The basemap toolkit requires python 2.4.""")
from distutils.core import setup

packages          = ['matplotlib.toolkits.basemap.data']
package_dirs       = {'':'lib'}
boundaryfiles = glob.glob("lib/matplotlib/toolkits/basemap/data/*_f.dat")
basemap_datafiles = [os.path.basename(bfile) for bfile in boundaryfiles]
package_data = {'matplotlib.toolkits.basemap.data':basemap_datafiles}
setup(
  name              = "basemap-data-fullres",
  version           = "0.9.7",
  description       = "full-resolution boundary data for basemap",
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
  package_data = package_data
  )
