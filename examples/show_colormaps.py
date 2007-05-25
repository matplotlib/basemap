from pylab import *
from matplotlib.toolkits.basemap import cm
a=outerproduct(arange(0,1,0.01),ones(10))
figure(figsize=(10,7))
subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
maps=[m for m in cm.datad.keys() if not m.endswith("_r")]
maps.sort()
l=len(maps)+1
i=1
for m in maps:
    subplot(1,l,i)
    axis("off")
    imshow(a,aspect='auto',cmap=cm.__dict__[m],origin="lower")
    title(m,rotation=90,fontsize=10)
    i=i+1
show()
