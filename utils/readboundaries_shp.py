import numpy as np
try:
    # try system version first
    from shapefile import Reader
except:
    # try bundled version as fallback
    from mpl_toolkits.basemap.shapefile import Reader

lsd = 5

def quantize(data,least_significant_digit):
    """

quantize data to improve compression. data is quantized using 
around(scale*data)/scale, where scale is 2**bits, and bits is determined 
from the least_significant_digit. For example, if 
least_significant_digit=1, bits will be 4.

This function is pure python.

    """
    precision = pow(10.,-least_significant_digit)
    exp = np.log10(precision)
    if exp < 0:
        exp = int(np.floor(exp))
    else:
        exp = int(np.ceil(exp))
    bits = np.ceil(np.log2(pow(10.,-exp)))
    scale = pow(2.,bits)
    return np.around(scale*data)/scale

def get_coast_polygons(resolution):
    polymeta = []; polybounds = []
    for level in [1,2,3,4]:
        filename = 'GSHHS_shp/%s/GSHHS_%s_L%s' % (resolution, resolution, level)
        #filename = 'WDBII_shp/%s/WDBII_border_%s_L%s' % (resolution, resolution, level)
        print filename
        shf = Reader(filename)
        fields = shf.fields
        try:
            shf.shapeRecords()
        except:
            continue
        for shprec in shf.shapeRecords():
            shp = shprec.shape; rec = shprec.record
            parts = shp.parts.tolist()
            if parts != [0]: 
                print 'multipart polygon'
                raise SystemExit
            verts = shp.points
            lons, lats = list(zip(*verts))
            north = max(lats); south = min(lats)
            attdict={}
            for r,key in zip(rec,fields[1:]):
                attdict[key[0]]=r
            area = attdict['area']
            id = attdict['id']
            polymeta.append([level,area,south,north,len(lons),id])
            b = np.empty((len(lons),2),np.float32)
            b[:,0] = lons; b[:,1] = lats
            if lsd is not None:
                b = quantize(b,lsd)
            polybounds.append(b)
    return polybounds, polymeta

def get_wdb_boundaries(resolution,level,rivers=False):
    polymeta = []; polybounds = []
    if rivers:
        filename = 'WDBII_shp/%s/WDBII_river_%s_L%02i' % (resolution, resolution, level)
    else:
        filename = 'WDBII_shp/%s/WDBII_border_%s_L%s' % (resolution, resolution, level)
    print filename
    shf = Reader(filename)
    fields = shf.fields
    for shprec in shf.shapeRecords():
        shp = shprec.shape; rec = shprec.record
        parts = shp.parts.tolist()
        if parts != [0]: 
            print 'multipart polygon'
            raise SystemExit
        verts = shp.points
        lons, lats = list(zip(*verts))
        north = max(lats); south = min(lats)
        attdict={}
        for r,key in zip(rec,fields[1:]):
            attdict[key[0]]=r
        area = -1
        id = attdict['id']
        polymeta.append([-1,-1,south,north,len(lons),id])
        b = np.empty((len(lons),2),np.float32)
        b[:,0] = lons; b[:,1] = lats
        if lsd is not None:
            b = quantize(b,lsd)
        polybounds.append(b)
    return polybounds, polymeta

# read in coastline data (only those polygons whose area > area_thresh).
for resolution in ['c','l','i','h','f']:
    poly, polymeta = get_coast_polygons(resolution)
    f = open('gshhs_'+resolution+'.dat','wb')
    f2 = open('gshhsmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        typ = pm[0]; area = pm[1]; south = pm[2]; north = pm[3]; npts = pm[4]
        id = pm[5]
        bstring = p.tostring()
        f.write(bstring)
        f2.write('%s %s %s %9.5f %9.5f %s %s %s\n' % (typ, area, npts, south,\
            north, offset, len(bstring),id))
        offset = offset + len(bstring)
    f.close()
    f2.close()
raise SystemExit

for resolution in ['c','l','i','h','f']:
    poly, polymeta = get_wdb_boundaries(resolution,1)
    f = open('countries_'+resolution+'.dat','wb')
    f2 = open('countriesmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        typ = pm[0]; area = pm[1]; south = pm[2]; north = pm[3]; npts = pm[4]
        id = pm[5]
        bstring = p.tostring()
        f.write(bstring)
        f2.write('%s %s %s %9.5f %9.5f %s %s %s\n' % (typ, area, npts, south,\
            north, offset, len(bstring),id))
        offset = offset + len(bstring)
    f.close()
    f2.close()

for resolution in ['c','l','i','h','f']:
    poly, polymeta = get_wdb_boundaries(resolution,2)
    f = open('states_'+resolution+'.dat','wb')
    f2 = open('statesmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        typ = pm[0]; area = pm[1]; south = pm[2]; north = pm[3]; npts = pm[4]
        id = pm[5]
        bstring = p.tostring()
        f.write(bstring)
        f2.write('%s %s %s %9.5f %9.5f %s %s %s\n' % (typ, area, npts, south,\
            north, offset, len(bstring),id))
        offset = offset + len(bstring)
    f.close()
    f2.close()

for resolution in ['c','l','i','h','f']:
    f = open('rivers_'+resolution+'.dat','wb')
    f2 = open('riversmeta_'+resolution+'.dat','w')
    for level in range(1,12):
        poly, polymeta = get_wdb_boundaries(resolution,level,rivers=True)
        offset = 0
        for p,pm in zip(poly,polymeta):
            typ = pm[0]; area = pm[1]; south = pm[2]; north = pm[3]; npts = pm[4]
            id = pm[5]
            bstring = p.tostring()
            f.write(bstring)
            f2.write('%s %s %s %9.5f %9.5f %s %s %s\n' % (typ, area, npts, south,\
                north, offset, len(bstring),id))
            offset = offset + len(bstring)
    f.close()
    f2.close()
