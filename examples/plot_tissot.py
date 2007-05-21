import pylab as p
from matplotlib.toolkits.basemap import Basemap as Basemap
from matplotlib.patches import Polygon

# Tissot's Indicatrix (http://en.wikipedia.org/wiki/Tissot's_Indicatrix). 
# These diagrams illustrate the distortion inherent in all map projections.
# In conformal projections, where angles are conserved around every location, 
# the Tissot's indicatrix are all circles, with varying sizes. In equal-area 
# projections, where area proportions between objects are conserved, the 
# Tissot's indicatrix have all unit area, although their shapes and 
# orientations vary with location.

# adapted from http://www.perrygeo.net/wordpress/?p=4

# create new figure
fig=p.figure()
m = Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80,
            projection='cyl')
shp_info = m.readshapefile('tissot','tissot',drawbounds=True)
ax = p.gca()
for nshape,seg in enumerate(m.tissot):
    poly = Polygon(seg,facecolor='green',zorder=10)
    ax.add_patch(poly)
# draw meridians and parallels.
m.drawparallels(p.arange(-90,91,30),labels=[1,0,0,0])
m.drawmeridians(p.arange(-180,180,60),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents()
p.title('Tissot Diagram - Cylindrical Equal Area')
print 'plot Cylindrical Equidistant Equal Area Tissot diagram ...'

# create new figure
fig=p.figure()
m = Basemap(llcrnrlon=-180,llcrnrlat=-70,urcrnrlon=180,urcrnrlat=70,
            projection='merc',lat_ts=20)
shp_info = m.readshapefile('tissot','tissot',drawbounds=True)
ax = p.gca()
for nshape,seg in enumerate(m.tissot):
    poly = Polygon(seg,facecolor='green',zorder=10)
    ax.add_patch(poly)
# draw meridians and parallels.
m.drawparallels(p.arange(-90,91,30),labels=[1,0,0,0])
m.drawmeridians(p.arange(-180,180,60),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents()
p.title('Tissot Diagram - Mercator Conformal')
print 'plot Mercator Conformal Tissot diagram ...'

# create new figure
fig=p.figure()
m = Basemap(lon_0=-60,lat_0=45,projection='ortho')
shp_info = m.readshapefile('tissot','tissot',drawbounds=False)
ax = p.gca()
for nshape,seg in enumerate(m.tissot):
    xx,yy = zip(*seg)
    if max(xx) < 1.e20 and max(yy) < 1.e20:
        poly = Polygon(seg,facecolor='green',zorder=10)
        ax.add_patch(poly)
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(p.arange(-90,91,30))
m.drawmeridians(p.arange(-180,180,30))
p.title('Tissot Diagram - Orthographic')
m.drawmapboundary()
p.gca().set_frame_on(True)
print 'plot Orthographic Tissot diagram ...'

# create new figure
fig=p.figure()
m = Basemap(lon_0=270,lat_0=90,boundinglat=10,projection='npstere')
shp_info = m.readshapefile('tissot','tissot',drawbounds=True)
ax = p.gca()
for nshape,seg in enumerate(m.tissot):
    poly = Polygon(seg,facecolor='green',zorder=10)
    ax.add_patch(poly)
# draw meridians and parallels.
m.drawparallels(p.arange(-90,91,30),labels=[1,0,0,0])
m.drawmeridians(p.arange(-180,180,30),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents()
p.title('Tissot Diagram - North Polar Stereographic Conformal')
print 'plot North Polar Stereographic Conformal Tissot diagram ...'

# create new figure
fig=p.figure()
m = Basemap(lon_0=270,lat_0=90,boundinglat=10,projection='nplaea')
shp_info = m.readshapefile('tissot','tissot',drawbounds=True)
ax = p.gca()
for nshape,seg in enumerate(m.tissot):
    poly = Polygon(seg,facecolor='green',zorder=10)
    ax.add_patch(poly)
# draw meridians and parallels.
m.drawparallels(p.arange(-90,91,30),labels=[1,0,0,0])
m.drawmeridians(p.arange(-180,180,30),labels=[0,0,0,1])
m.drawcoastlines()
m.fillcontinents()
p.title('Tissot Diagram - North Polar Lambert Azimuthal Equal Area')
print 'plot North Polar Lambert Azimuthal Equal Area Tissot diagram ...'
p.show()
