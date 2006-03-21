# make plots of etopo bathymetry/topography data on
# various map projections, drawing coastlines, state and
# country boundaries, filling continents and drawing
# parallels/meridians

# illustrates special-case polar-centric projections.

from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topodatin = array(load('etopo20data.gz'),'d')
lonsin = array(load('etopo20lons.gz'),'d')
latsin = array(load('etopo20lats.gz'),'d')

# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topodatin,lonsin,start=False)
lats = latsin

print 'min/max etopo20 data:'
print min(ravel(topoin)),max(ravel(topoin))

boundinglat = 20. 
projs = ['laea','stere','aeqd']
projnames = ['Lambert Azimuthal Equal Area','Stereographic','Azimuthal Equidistant']
for proj,projname in zip(projs,projnames):
    # setup stereographic map projection (Southern Hemisphere).
    # centered on Australia
    lon_0 = 130.
    m = Basemap(boundinglat=-boundinglat,lon_0=lon_0,\
        	resolution='c',area_thresh=10000.,projection='sp'+proj)
    # transform to nx x ny regularly spaced native projection grid
    nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
    topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
    # setup figure with same aspect ratio as map.
    fig = figure(figsize=(10,6))
    ax = fig.add_subplot(121)
    # plot image over map.
    im = m.imshow(topodat,cm.jet)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    # draw parallels and meridians (labelling is 
    # not implemented for orthographic).
    parallels = arange(-80.,90,20.)
    m.drawparallels(parallels)
    meridians = arange(0.,360.,60.)
    m.drawmeridians(meridians)
    title('South Polar '+projname,y=1.075)

    # setup of basemap ('ortho' = orthographic projection)
    m = Basemap(projection='ortho',
        	resolution='c',area_thresh=10000.,lat_0=-90,lon_0=lon_0-180.)
    ax = fig.add_subplot(122)
    x,y = m(*meshgrid(lonsin,latsin))
    cs = m.contourf(x,y,topodatin,20,cmap=cm.jet)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    # draw parallels and meridians (labelling is 
    # not implemented for orthographic).
    parallels = arange(-80.,90,20.)
    m.drawparallels(parallels)
    meridians = arange(0.,360.,60.)
    m.drawmeridians(meridians)
    # draw boundary around map region.
    m.drawmapboundary()
    title('South Polar Orthographic',y=1.075)
    show()

    # setup stereographic map projection (Northern Hemisphere).
    # centered on US
    lon_0 = -90.
    m = Basemap(boundinglat=boundinglat,lon_0=lon_0,\
        	resolution='c',area_thresh=10000.,projection='np'+proj)
    # transform to nx x ny regularly spaced native projection grid
    nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
    topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
    # setup figure with same aspect ratio as map.
    fig = figure(figsize=(10,6))
    ax = fig.add_subplot(121)
    # plot image over map.
    im = m.imshow(topodat,cm.jet)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    # draw parallels and meridians (labelling is 
    # not implemented for orthographic).
    parallels = arange(-80.,90,20.)
    m.drawparallels(parallels)
    meridians = arange(0.,360.,60.)
    m.drawmeridians(meridians)
    title('North Polar '+projname,y=1.075)

    # setup of basemap ('ortho' = orthographic projection)
    m = Basemap(projection='ortho',
        	resolution='c',area_thresh=10000.,lat_0=90,lon_0=lon_0)
    ax = fig.add_subplot(122)
    x,y = m(*meshgrid(lonsin,latsin))
    cs = m.contourf(x,y,topodatin,20,cmap=cm.jet)
    # draw coastlines and political boundaries.
    m.drawcoastlines()
    # draw parallels and meridians (labelling is 
    # not implemented for orthographic).
    parallels = arange(-80.,90,20.)
    m.drawparallels(parallels)
    meridians = arange(0.,360.,60.)
    m.drawmeridians(meridians)
    # draw boundary around map region.
    m.drawmapboundary()
    title('North Polar Orthographic',y=1.075)
    show()
