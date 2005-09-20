from distutils.core import setup
import os

datadir = os.environ.get('BASEMAP_DATA_PATH')
if not datadir:
    datadir ='share/basemap'

setup(
  version           = "0.7.1",
  description       = "boundary data for basemap matplotlib toolkit",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  download_url      = "http://sourceforge.net/projects/matplotlib",
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  license           = "OSI Approved",
  name              = "basemap-data",
  data_files        = [(datadir,['data/countries_c.txt','data/states_c.txt','data/countries_l.txt','data/states_l.txt','data/gshhs_c.txt','data/gshhs_l.txt','data/countries_i.txt','data/states_i.txt','data/gshhs_i.txt','data/countries_h.txt','data/states_h.txt','data/gshhs_h.txt','data/rivers_c.txt','data/rivers_l.txt','data/rivers_i.txt','data/rivers_h.txt'])])
