import os
import numpy as np
from shapefile import Reader

lsd = 5

UTILS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(UTILS_DIR, '..', 'lib', 'mpl_toolkits',
                          'basemap', 'data')

# Folder where GSHHG shapefiles were extracted. Change if needed
GSHHS_DIR = UTILS_DIR

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

def interpolate_long_segments(coords, resolution):
    lookup_thresh = {'c': 0.5, 'l':0.3, 'i':0.2, 'h':0.1, 'f':0.05}
    thresh = lookup_thresh[resolution]
    spacing = thresh / 5.0

    lons, lats = coords.T
    dist = np.hypot(np.diff(lons), np.diff(lats))

    if np.all(dist <= thresh):
        return coords

    out_lon, out_lat = [], []
    for i in np.arange(len(dist)):
        if dist[i] <= thresh:
            out_lon.append(lons[i])
            out_lat.append(lats[i])
        else:
            x = [0, dist[i]]
            new_x = np.arange(0, dist[i], spacing)
            out_lon.extend(np.interp(new_x, x, lons[i:i+2]))
            out_lat.extend(np.interp(new_x, x, lats[i:i+2]))

    out_lon.append(lons[-1])
    out_lat.append(lats[-1])
    return np.column_stack([out_lon, out_lat]).astype(coords.dtype)

def get_coast_polygons(resolution):
    polymeta = []; polybounds = []
    for level in [1,2,3,5]:
        filename = os.path.join(GSHHS_DIR, 'GSHHS_shp/', resolution,
                                'GSHHS_{}_L{}'.format(resolution, level))
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

        # Manual fix for incorrect Antarctica polygons at full resolution
        # This issue is only present in the shapefile version and may be fixed
        # in future versions of GSHHS!
        if resolution == 'f' and level == 5:
            i = [item[-1] for item in polymeta].index('4-E')
            coords = polybounds[i][2:-1, :]
            coords = np.vstack([coords,
                                [180.0, -90.0],
                                [0.0, -90.0]]).astype(np.float32)
            polybounds[i] = coords
            polymeta[i][-2] = len(coords)

            j = [item[-1] for item in polymeta].index('4-W')
            coords = polybounds[j][3:, :]
            np.savetxt('coordinates.txt', coords)
            coords = np.vstack([coords,
                                [0.0, coords[-1][1]],
                                [0.0, -90.0],
                                [-180.0, -90.0],
                                coords[0]]).astype(np.float32)

            polybounds[j] = coords
            polymeta[j][-2] = len(coords)

    return polybounds, polymeta

def get_wdb_boundaries(resolution,level,rivers=False):
    polymeta = []; polybounds = []
    if rivers:
        filename = os.path.join(GSHHS_DIR, 'WDBII_shp', resolution,
                            'WDBII_river_{}_L{:02}'.format(resolution, level))
    else:
        filename = os.path.join(GSHHS_DIR, 'WDBII_shp', resolution,
                            'WDBII_border_{}_L{}'.format(resolution, level))
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
        # Detect degenerate lines that are actually points...
        if len(verts) == 2 and np.allclose(verts[0], verts[1]):
            print 'Skipping degenerate line...'
            continue

        lons, lats = list(zip(*verts))
        north = max(lats); south = min(lats)
        attdict={}
        for r,key in zip(rec,fields[1:]):
            attdict[key[0]]=r
        area = -1
        poly_id = attdict['id']
        b = np.empty((len(lons),2),np.float32)
        b[:,0] = lons; b[:,1] = lats

        if not rivers:
            b = interpolate_long_segments(b, resolution)

        if lsd is not None:
            b = quantize(b,lsd)

        polymeta.append([-1,-1,south,north,len(b),poly_id])
        polybounds.append(b)

    return polybounds, polymeta

# read in coastline data (only those polygons whose area > area_thresh).
for resolution in ['c','l','i','h','f']:
    poly, polymeta = get_coast_polygons(resolution)
    f = open(os.path.join(OUTPUT_DIR, 'gshhs_'+resolution+'.dat'), 'wb')
    f2 = open(os.path.join(OUTPUT_DIR, 'gshhsmeta_'+resolution+'.dat'), 'w')
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
    poly, polymeta = get_wdb_boundaries(resolution,1)
    f = open(os.path.join(OUTPUT_DIR, 'countries_'+resolution+'.dat'), 'wb')
    f2 = open(os.path.join(OUTPUT_DIR, 'countriesmeta_'+resolution+'.dat'), 'w')
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
    f = open(os.path.join(OUTPUT_DIR, 'states_'+resolution+'.dat'), 'wb')
    f2 = open(os.path.join(OUTPUT_DIR, 'statesmeta_'+resolution+'.dat'), 'w')
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
    f = open(os.path.join(OUTPUT_DIR, 'rivers_'+resolution+'.dat'), 'wb')
    f2 = open(os.path.join(OUTPUT_DIR, 'riversmeta_'+resolution+'.dat'), 'w')
    # Levels above 5 are intermittent rivers and irrigation canals.
    # They haven't been included in the past, as far as I can tell, so I'm
    # not including them here...
    offset = 0
    for level in range(1, 6):
        poly, polymeta = get_wdb_boundaries(resolution,level,rivers=True)
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
