from distutils.core import setup
import os

datadir = os.environ.get('BASEMAP_DATA_PATH')
if not datadir:
    datadir ='share/basemap'

setup(
  version           = "0.7",
  description       = "boundary data for basemap matplotlib toolkit",
  url               = "http://matplotlib.sourceforge.net/toolkits.html",
  download_url      = "http://sourceforge.net/projects/matplotlib",
  author            = "Jeff Whitaker",
  author_email      = "jeffrey.s.whitaker@noaa.gov",
  license           = "OSI Approved",
  name              = "basemap-data",
  data_files        = [(datadir,['countries_c.txt','states_c.txt','countries_l.txt','states_l.txt','gshhs_c.txt','gshhs_l.txt','countries_i.txt','states_i.txt','gshhs_i.txt','countries_h.txt','states_h.txt','gshhs_h.txt','rivers_c.txt','rivers_l.txt','rivers_i.txt','rivers_h.txt'])])
