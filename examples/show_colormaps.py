import numpy, pylab
from mpl_toolkits.basemap import cm
a=numpy.outer(numpy.arange(0,1,0.01),numpy.ones(10))
pylab.figure(figsize=(10,7))
pylab.subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
maps=[m for m in cm.datad.keys() if not m.endswith("_r")]
maps.sort()
l=len(maps)+1
i=1
for m in maps:
    pylab.subplot(1,l,i)
    pylab.axis("off")
    pylab.imshow(a,aspect='auto',cmap=cm.__dict__[m],origin="lower")
    pylab.title(m,rotation=90,fontsize=10)
    i=i+1
pylab.show()
