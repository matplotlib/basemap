import sys, glob, os
if 'setuptools' in sys.modules:
# Are we running with setuptools?
# if so, need to specify all the packages in heirarchy
    additional_params = {'namespace_packages' : ['matplotlib.toolkits']}    
    packages.extend(['matplotlib', 'matplotlib.toolkits'])
    setup = setuptools.setup
else:
    additional_params = {}
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
