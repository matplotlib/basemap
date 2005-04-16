import matplotlib
from matplotlib.toolkits.basemap import Basemap
from pylab import *

m = Basemap(-11.,51.,-5.,56.,
            resolution='c',area_thresh=1000.,projection='cyl')
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
ax = gca() # get current axis instance
ax.update_datalim(((m.llcrnrx, m.llcrnry),(m.urcrnrx,m.urcrnry)))
ax.set_xlim((m.llcrnrx, m.urcrnrx))
ax.set_ylim((m.llcrnry, m.urcrnry))
m.drawcoastlines()
m.fillcontinents()
# draw parallels
circles = arange(50,60,1).tolist()
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
meridians = arange(-12,0,1)
m.drawmeridians(meridians,labels=[1,1,1,1])
title("Crude Res Coastlines ('c')",y=1.075)
show()

m = Basemap(-11.,51.,-5.,56.,
            resolution='l',area_thresh=1000.,projection='cyl')
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
ax = gca() # get current axis instance
ax.update_datalim(((m.llcrnrx, m.llcrnry),(m.urcrnrx,m.urcrnry)))
ax.set_xlim((m.llcrnrx, m.urcrnrx))
ax.set_ylim((m.llcrnry, m.urcrnry))
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title("Low Res Coastlines ('l')",y=1.075)
show()

m = Basemap(-11.,51.,-5.,56.,
            resolution='i',area_thresh=1000.,projection='cyl')
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
ax = gca() # get current axis instance
ax.update_datalim(((m.llcrnrx, m.llcrnry),(m.urcrnrx,m.urcrnry)))
ax.set_xlim((m.llcrnrx, m.urcrnrx))
ax.set_ylim((m.llcrnry, m.urcrnry))
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title("Intermediate Res Coastlines ('i')",y=1.075)
show()
