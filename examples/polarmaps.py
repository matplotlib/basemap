from __future__ import print_function
# make plots of etopo bathymetry/topography data on
# various map projections, drawing coastlines, state and
# country boundaries, filling continents and drawing
# parallels/meridians

# illustrates special-case polar-centric projections.

from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
etopo = np.loadtxt('etopo20data.gz')
lons = np.loadtxt('etopo20lons.gz')
lats = np.loadtxt('etopo20lats.gz')

print('min/max etopo20 data:')
print(etopo.min(),etopo.max())

# these are the 4 polar projections
projs = ['laea','stere','aeqd','ortho'] # short names
# long names
projnames = ['Lambert Azimuthal Equal Area','Stereographic','Azimuthal Equidistant','Orthographic']
# loop over hemispheres, make a 4-panel plot for each hemisphere
# showing all four polar projections.
for hem in ['North','South']:
    if hem == 'South':
        lon_0 = -130.
        lon_0_ortho = lon_0 - 180.
        lat_0 = -90.
        #  Lambert Azimuth bounding lat must not extend into opposite hem.
        bounding_lat = -0.01
    elif hem == 'North':
        lon_0 = 130.
        lon_0_ortho = lon_0
        lat_0 = 90.
        #  Lambert Azimuth bounding lat must not extend into opposite hem.
        bounding_lat = 0.01
    # loop over projections, one for each panel of the figure.
    fig = plt.figure(figsize=(8,8))
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
                       resolution='c',area_thresh=10000.,projection=projection,round=True)
        # compute native map projection coordinates for lat/lon grid.
        x,y = m(*np.meshgrid(lons,lats))
        ax = fig.add_subplot(2,2,npanel)
        # make filled contour plot.
        cs = m.contourf(x,y,etopo,np.linspace(-7500,4500,41),cmap=cm.GMT_haxby)
        # draw coastlines.
        m.drawcoastlines()
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,90,20.))
        #labels = [l,r,t,b]
        m.drawmeridians(np.arange(0.,340.,30.),labels=[1,1,1,1],fontsize=7)
        # draw boundary around map region.
        m.drawmapboundary()
        # draw title.
        plt.title(hem+' Polar '+projname,y=1.05,fontsize=12)
        print('plotting '+hem+' Polar '+projname+' basemap ...')
plt.show()
