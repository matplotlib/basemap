from distutils.core import setup, Extension
import sys, glob

deps = glob.glob('src/*.c')

extensions = [Extension("proj4",deps,include_dirs = ['src'],)]

datadir ='share/basemap-py'+repr(sys.version_info[0])+repr(sys.version_info[1])

setup(
  name              = "basemap",
  version           = "0.5.3cvs",
  description       = "Plot data on map projections with matplotlib",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  data_files        = [(datadir,['data/countries_c.txt','data/states_c.txt','data/countries_l.txt','data/states_l.txt','data/gshhs_c.txt','data/gshhs_l.txt','data/countries_i.txt','data/states_i.txt','data/gshhs_i.txt'])],
  packages          = ['matplotlib/toolkits','matplotlib/toolkits/basemap'],
  package_dir       = {'':'lib'},
  ext_modules       = extensions)
