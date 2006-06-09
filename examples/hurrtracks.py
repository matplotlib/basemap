"""
draw Atlantic Hurricane Tracks for storms that reached Cat 4 or 5.
part of the track for which storm is cat 4 or 5 is shown red.
ESRI shapefile data from http://www.nationalatlas.gov/atlasftp.html
"""
import pylab as p
from matplotlib.toolkits.basemap import Basemap as Basemap
# Lambert Conformal Conic map.
m = Basemap(llcrnrlon=-100.,llcrnrlat=0.,urcrnrlon=-20.,urcrnrlat=57.,
            projection='lcc',lat_1=20.,lat_2=40.,lon_0=-60.,
            resolution ='l',area_thresh=1000.)
# create figure, add axes.
fig=p.figure()
fig.add_axes([0.1,0.1,0.8,0.8],axisbg='#99ffff')
# read shapefile.
shp_info = m.readshapefile('huralll020','hurrtracks',drawbounds=False)
print shp_info
# find names of storms that reached Cat 4.
names = []
for shapedict in m.hurrtracks_info:
    cat = shapedict['CATEGORY']
    name = shapedict['NAME']
    if cat in ['H4','H5'] and name not in names:
        # only use named storms.
        if name != 'NOT NAMED':  names.append(name)
print names
print len(names)
# plot tracks of those storms.
for shapedict,shape in zip(m.hurrtracks_info,m.hurrtracks):
    name = shapedict['NAME']
    cat = shapedict['CATEGORY']
    if name in names:
        xx,yy = zip(*shape)
        # show part of track where storm > Cat 4 as thick red.
        if cat in ['H4','H5']: 
            p.plot(xx,yy,linewidth=1.5,color='r')
        elif cat in ['H1','H2','H3']:
            p.plot(xx,yy,color='k')
# draw coastlines, meridians and parallels.
m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='#cc9966')
m.drawparallels(p.arange(10,70,20),labels=[1,1,0,0])
m.drawmeridians(p.arange(-100,0,20),labels=[0,0,0,1])
p.title('Atlantic Hurricane Tracks (Storms Reaching Category 4, 1851-2004)')
p.show()
