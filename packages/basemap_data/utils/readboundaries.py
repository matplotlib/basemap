import sys
import numpy as np

lsd = 3

def quantize(data,least_significant_digit):
    """
    Quantize data to improve compression. data is quantized using
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

def get_coast_polygons(coastfile):
    polymeta = []; polybounds = []
    lats = []; lons = []
    for line in open(coastfile):
        if line.startswith('#'):
            continue
        linesplit = line.strip().split()
        if line.startswith('>'):
            area, west, east, south, north = map(float, linesplit[5:10])
            poly_id = linesplit[-1]
            level = linesplit[3]
            polymeta.append([level,area,south,north,poly_id])
            if lons:
                #lons.append(lons[0]); lats.append(lats[0])
                b = np.empty((len(lons),2),np.float32)
                b[:,0] = lons; b[:,1] = lats
                if lsd is not None:
                    b = quantize(b,lsd)
                polybounds.append(b)
            lats = []; lons = []
            continue
        lon = float(linesplit[0])
        lat = float(linesplit[1])
        lons.append(lon); lats.append(lat)
    #lons.append(lons[0]); lats.append(lats[0])
    b = np.empty((len(lons),2),np.float32)
    b[:,0] = lons; b[:,1] = lats
    if lsd is not None:
        b = quantize(b,lsd)
    polybounds.append(b)
    polymeta2 = []
    for meta,bounds in zip(polymeta,polybounds):
        npts = bounds.shape[0]
        polymeta2.append(meta[:-1] + [npts] + [meta[-1]])
    return polybounds, polymeta2

def get_boundary_lines(bdatfile, resolution):
    lons = []; lats = []; polybounds = []
    for line in open(bdatfile):
        if line.startswith('#'): continue
        linesplit = line.split()
        if line.startswith('>'):
           if lons:
               b = np.empty((len(lons),2),np.float32)
               b[:,0] = lons; b[:,1] = lats
               b = interpolate_long_segments(b, resolution)
               if lsd is not None:
                   b = quantize(b,lsd)
               polybounds.append(b)
           lons = []; lats = []
           continue
        lon, lat = [float(val) for val in linesplit]
        lats.append(lat); lons.append(lon)
    b = np.empty((len(lons),2),np.float32)
    b[:,0] = lons; b[:,1] = lats
    b = interpolate_long_segments(b, resolution)
    if lsd is not None:
        b = quantize(b,lsd)
    polybounds.append(b)
    polymeta = []
    polybounds2 = []
    for bounds in polybounds:
        npts = bounds.shape[0]
        if npts == 2 and\
           bounds[0,0] == bounds[1,0] and\
           bounds[0,1] == bounds[1,1]: continue
        polybounds2.append(bounds)
        south = bounds[:,1].min()
        north = bounds[:,1].max()
        polymeta.append((south,north,npts))
    return polybounds2, polymeta

# read in coastline data (only those polygons whose area > area_thresh).
for resolution in ['c','l','i','h','f']:
    coastlons = []; coastlats = []; coastsegind = []; coastsegtype = []
    coastfile = 'gshhs_'+resolution+'.txt'
    countryfile = 'countries_'+resolution+'.txt'
    statefile = 'states_'+resolution+'.txt'
    riverfile = 'rivers_'+resolution+'.txt'

    poly, polymeta = get_coast_polygons(coastfile)
    f = open('../lib/mpl_toolkits/basemap/data/gshhs_'+resolution+'.dat','wb')
    f2 = open('../lib/mpl_toolkits/basemap/data/gshhsmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        typ = pm[0]; area = pm[1]; south = pm[2]; north = pm[3]; npts = pm[4]
        poly_id = pm[5]
        f2.write('%s %s %s %9.5f %9.5f %s %s %s\n' % (typ, area, npts, south, north, offset, len(bstring), poly_id))
        offset = offset + len(bstring)
    f.close()
    f2.close()

    poly, polymeta = get_boundary_lines(countryfile, resolution)
    f = open('../lib/mpl_toolkits/basemap/data/countries_'+resolution+'.dat','wb')
    f2 = open('../lib/mpl_toolkits/basemap/data/countriesmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        south,north,npts = pm[:]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (-1,-1,npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()

    poly, polymeta = get_boundary_lines(statefile, resolution)
    f = open('../lib/mpl_toolkits/basemap/data/states_'+resolution+'.dat','wb')
    f2 = open('../lib/mpl_toolkits/basemap/data/statesmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        south,north,npts = pm[:]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (-1,-1,npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()

    poly, polymeta = get_boundary_lines(riverfile, resolution)
    f = open('../lib/mpl_toolkits/basemap/data/rivers_'+resolution+'.dat','wb')
    f2 = open('../lib/mpl_toolkits/basemap/data/riversmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        south,north,npts = pm[:]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (-1,-1,npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()
