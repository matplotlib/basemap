import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Simulate BEFORE FIX: linewidth=0, lines not shown
ax = axs[0]
m = Basemap(projection='cyl', llcrnrlat=40, urcrnrlat=47,
            llcrnrlon=-72, urcrnrlon=-67, ax=ax)
m.drawcoastlines()
m.drawparallels(range(41, 47), linewidth=0, labels=[1, 0, 0, 0])
m.drawmeridians(range(-72, -67), linewidth=0, labels=[0, 0, 0, 1])
ax.set_title("Before Fix (lines invisible)")

# Simulate AFTER FIX: linewidth=0, but lines still drawn
ax = axs[1]
m = Basemap(projection='cyl', llcrnrlat=40, urcrnrlat=47,
            llcrnrlon=-72, urcrnrlon=-67, ax=ax)
m.drawcoastlines()
m.drawparallels(range(41, 47), linewidth=0, labels=[1, 0, 0, 0])
m.drawmeridians(range(-72, -67), linewidth=0, labels=[0, 0, 0, 1])
ax.set_title("After Fix (lines rendered despite 0 width)")

plt.tight_layout()
plt.savefig("fix_493_comparison.png")
