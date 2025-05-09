import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

# Common settings
llcrnrlon, llcrnrlat = -72, 41
urcrnrlon, urcrnrlat = -67, 47
parallels = np.arange(41., 47., 1.)

# BEFORE FIX (simulate broken dashes)
ax1 = axes[0]
ax1.set_title("Before Fix (lines invisible)")
m1 = Basemap(projection='cyl', resolution='l',
             llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
             urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, ax=ax1)
m1.drawcoastlines()
for lat in parallels:
    x, y = m1(np.linspace(llcrnrlon, urcrnrlon, 100), lat * np.ones(100))
    l = plt.Line2D(x, y, linewidth=1.0)
    l.set_dashes([1e-10, 1e-10])  # Buggy dash
    ax1.add_line(l)

# AFTER FIX (skip broken dashes)
ax2 = axes[1]
ax2.set_title("After Fix (lines visible)")
m2 = Basemap(projection='cyl', resolution='l',
             llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
             urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat, ax=ax2)
m2.drawcoastlines()
for lat in parallels:
    x, y = m2(np.linspace(llcrnrlon, urcrnrlon, 100), lat * np.ones(100))
    l = plt.Line2D(x, y, linewidth=1.0)
    l.set_dashes([1, 1])  # Fixed dash
    ax2.add_line(l)

plt.tight_layout()
plt.savefig("dash_fix_comparison.pdf")
print("Saved as dash_fix_comparison.pdf")
