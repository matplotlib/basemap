# make plots of etopo bathymetry/topography data on
# various map projections, drawing coastlines, state and
# country boundaries, filling continents and drawing
# parallels/meridians

# illustrates special-case polar-centric projections.

from matplotlib.toolkits.basemap import Basemap
from pylab import title, colorbar, show, axes, cm, load, arange, \
                  figure, ravel, meshgrid

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
etopo = load('etopo20data.gz')
lons = load('etopo20lons.gz')
lats = load('etopo20lats.gz')

print 'min/max etopo20 data:'
print min(ravel(etopo)),max(ravel(etopo))

# these are the 4 polar projections
projs = ['laea','stere','aeqd','ortho'] # short names
# long names
projnames = ['Lambert Azimuthal Equal Area','Stereographic','Azimuthal Equidistant','Orthographic']
# loop over hemispheres, make a 4-panel plot for each hemisphere
# showing all four polar projections.
for hem in ['North','South']:
    if hem == 'South':
        lon_0 = 130.
        lon_0_ortho = lon_0 - 180.
        lat_0 = -90.
        bounding_lat = -20.
    elif hem == 'North':
        lon_0 = -90.
        lon_0_ortho = lon_0
        lat_0 = 90.
        bounding_lat = 20.
    # loop over projections, one for each panel of the figure.
    fig = figure(figsize=(8,8))
    npanel = 0
    for proj,projname in zip(projs,projnames):
        npanel = npanel + 1
        if hem == 'South':
            projection = 'sp'+proj
        elif hem == 'North':
            projection = 'np'+proj
        # setup map projection
        # centered on Australia (for SH) or US (for NH).
        if proj == 'ortho':
           m = Basemap(projection='ortho',
                       resolution='c',area_thresh=10000.,lat_0=lat_0,lon_0=lon_0_ortho)
        else:
           m = Basemap(boundinglat=bounding_lat,lon_0=lon_0,\
                       resolution='c',area_thresh=10000.,projection=projection)
        # compute native map projection coordinates for lat/lon grid.
        x,y = m(*meshgrid(lons,lats))
        ax = fig.add_subplot(2,2,npanel)
        # make filled contour plot.
        cs = m.contourf(x,y,etopo,20,cmap=cm.jet)
        # draw coastlines.
        m.drawcoastlines()
        # draw parallels and meridians.
        m.drawparallels(arange(-80.,90,20.))
        m.drawmeridians(arange(0.,360.,60.))
        # draw boundary around map region.
        m.drawmapboundary()
        # draw title.
        title(hem+' Polar '+projname,y=1.05,fontsize=12)
        print 'plotting '+hem+' Polar '+projname+' basemap ...'
show()
