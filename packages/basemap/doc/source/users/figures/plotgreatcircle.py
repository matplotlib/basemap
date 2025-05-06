from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

def better_great_circle(m, lon1, lat1, lon2, lat2, arc_height_ratio=-1.5, **kwargs):
    # Interpolate linearly
    t = np.linspace(0, 1, 300)
    lons = lon1 + (lon2 - lon1) * t
    lats = lat1 + (lat2 - lat1) * t

    # Determine arc direction: raise toward equator
    mid_lat = (lat1 + lat2) / 2
    arc_direction = -1 if mid_lat > 0 else 1  # bend toward equator

    # Create realistic arc using cosine curve
    arc_height = arc_height_ratio * abs(lat2 - lat1 + 1e-6)  # prevent divide by 0
    bulge = np.cos(np.pi * (t - 0.5))  # peak at center
    lats += arc_direction * arc_height * bulge

    # Clip to Mercator-safe bounds
    lats = np.clip(lats, -85, 85)

    # Project and plot
    x, y = m(lons, lats)
    mask = np.isfinite(x) & np.isfinite(y)
    m.plot(x[mask], y[mask], **kwargs)


# create new figure, axes instances.
fig=plt.figure(figsize=(14, 10))
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup mercator map projection.

# Change it to full Map view so we can see other great circle route
m = Basemap(projection='merc',
            llcrnrlon=-180, urcrnrlon=180,
            llcrnrlat=-60, urcrnrlat=80,
            lat_ts=20, resolution='c')

m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(10, 90, 20),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])

# nylat, nylon are lat/lon of New York
nylat = 40.78; nylon = -73.98
# lonlat, lonlon are lat/lon of London.
lonlat = 51.53; lonlon = 0.08
# draw great circle route between NY and London
m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')

# # Anchorage → Moscow with gap
# m.drawgreatcircle(-149.9003, 61.2181, 37.6173, 55.7558, color='red', linewidth=2)
#
# # Nome → Helsinki with gap
# m.drawgreatcircle(-165.4064, 64.5011, 24.9354, 60.1695, color='blue', linewidth=2)
#
# # Fairbanks → Murmansk with gap
# m.drawgreatcircle(-147.7164, 64.8378, 33.0846, 68.9585,  color='green', linewidth=2)

# Anchorage → Moscow
better_great_circle(m, -149.9003, 61.2181, 37.6173, 55.7558, color='red', linewidth=2)

# Nome → Helsinki
better_great_circle(m, -165.4064, 64.5011, 24.9354, 60.1695, color='blue', linewidth=2)

# Fairbanks → Murmansk
better_great_circle(m, -147.7164, 64.8378, 33.0846, 68.9585, color='green', linewidth=2)

# Great circle: San Francisco → Tokyo
# m.drawgreatcircle(-122.4194, 37.7749, 139.6917, 35.6895, linewidth=2, color='b')

# Great circle: Cape Town → Beijing
# m.drawgreatcircle(18.4241, -33.9249, 116.4074, 39.9042, linewidth=2, color='r')

# Great circle: Rio → Sydney
# m.drawgreatcircle(-43.1729, -22.9068, 151.2093, -33.8688, linewidth=2, color='g')

ax.set_title("Fixed Great Circles on Mercator Projection")
plt.show()














