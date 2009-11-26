from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# 'cubed sphere'
# inscribe the sphere in a cube, then separately project each cube
# face with gnomonic projection.
# http://www.progonos.com/furuti/MapProj/Normal/ProjPoly/Foldout/Cube/cube.html
# suitable for cutting and folding.

# choose figure size to match aspect ratio of map.
fig = plt.figure(figsize=(10,7.5))
fig.subplots_adjust(bottom=0, left=0, right=1, top=1, wspace=0, hspace=0)
rsphere = 6370997.
width = 2.*rsphere; height=width
npanel=0
for lat_0 in [90,0,-90]:
    for ncol in range(0,4):
        npanel = npanel + 1
        if lat_0 == 0 or ncol == 1:
            ax=fig.add_subplot(3,4,npanel)
            ax.set_frame_on(False)
            lon_0=225 + 90*(ncol+1) - 45
            # use fix_aspect=False, so white space won't appear between
            # faces of cube when window is resized.
            m = Basemap(width=width,height=height,resolution=None,\
                        projection='gnom',lon_0=lon_0,lat_0=lat_0,\
                        rsphere=rsphere,fix_aspect=False)
            m.bluemarble(scale=0.5)
            m.drawparallels(np.arange(-90,91,10),color='0.5')
            m.drawmeridians(np.arange(5,365,10),color='0.5')
            #m.drawlsmask(ocean_color='aqua',land_color='coral')
            #m.drawparallels(np.arange(-90,91,10))
            #m.drawmeridians(np.arange(5,365,10))
fig.text(0.625,0.75,\
        'World Map on a Cube\n Gnomonic Projection',\
        fontsize=14)
plt.show()
