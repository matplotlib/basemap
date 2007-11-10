import numpy, sys

lsd = 3

def quantize(data,least_significant_digit):
    """

quantize data to improve compression. data is quantized using 
around(scale*data)/scale, where scale is 2**bits, and bits is determined 
from the least_significant_digit. For example, if 
least_significant_digit=1, bits will be 4.

This function is pure python.

    """
    precision = pow(10.,-least_significant_digit)
    exp = numpy.log10(precision)
    if exp < 0:
        exp = int(numpy.floor(exp))
    else:
        exp = int(numpy.ceil(exp))
    bits = numpy.ceil(numpy.log2(pow(10.,-exp)))
    scale = pow(2.,bits)
    return numpy.around(scale*data)/scale

def get_coast_polygons(coastfile,area_thresh=-1.,latmin=-91.,latmax=91.):
    polymeta = []; polybounds = []
    lats = []; lons = []
    for line in open(coastfile):
        linesplit = line.split()
        if line.startswith('P'):
            area = float(linesplit[5])
            west,east,south,north = float(linesplit[6]),float(linesplit[7]),float(linesplit[8]),float(linesplit[9])
            typ = int(linesplit[3])
            useit = latmax>=south and latmin<=north and area>area_thresh
            if useit:
                polymeta.append((typ,area,south,north))
                if lons:
                    #lons.append(lons[0]); lats.append(lats[0])
                    b = numpy.empty((len(lons),2),numpy.float32)
                    b[:,0] = lons; b[:,1] = lats
                    if lsd is not None:
                        b = quantize(b,lsd)
                    polybounds.append(b)
                lats = []; lons = []
            continue
        if useit:
            lon, lat = [float(val) for val in linesplit]
            lons.append(lon); lats.append(lat)
    #lons.append(lons[0]); lats.append(lats[0])
    b = numpy.empty((len(lons),2),numpy.float32)
    b[:,0] = lons; b[:,1] = lats
    if lsd is not None:
        b = quantize(b,lsd)
    polybounds.append(b)
    polymeta2 = []
    for meta,bounds in zip(polymeta,polybounds):
        npts = bounds.shape[0]
        polymeta2.append(meta+(npts,))
    return polybounds, polymeta2

def get_boundary_lines(bdatfile,latmin=-91.,latmax=91.):
    polymeta = []; polybounds = []
    lons = []; lats = []
    for line in open(bdatfile):
        linesplit = line.split()
        if line.startswith('>'):
            west,east,south,north = float(linesplit[7]),float(linesplit[8]),float(linesplit[9]),float(linesplit[10])
            useit = latmax>=south and latmin<=north
            if useit:
                polymeta.append((south,north))
                if lons:
                    b = numpy.empty((len(lons),2),numpy.float32)
                    b[:,0] = lons; b[:,1] = lats
                    if lsd is not None:
                        b = quantize(b,lsd)
                    polybounds.append(b)
                lons = []; lats = []
            continue
        if useit:
            lon, lat = [float(val) for val in linesplit]
            lats.append(lat); lons.append(lon)
    b = numpy.empty((len(lons),2),numpy.float32)
    b[:,0] = lons; b[:,1] = lats
    if lsd is not None:
        b = quantize(b,lsd)
    polybounds.append(b)
    polymeta2 = []
    polybounds2 = []
    for meta,bounds in zip(polymeta,polybounds):
        npts = bounds.shape[0]
        if npts == 2 and\
           bounds[0,0] == bounds[1,0] and\
           bounds[0,1] == bounds[1,1]: continue
        polybounds2.append(bounds)
        polymeta2.append(meta+(npts,))
    return polybounds2, polymeta2

# read in coastline data (only those polygons whose area > area_thresh).
#resolution = sys.argv[1]
for resolution in ['c','l','i','h']:
    coastlons = []; coastlats = []; coastsegind = []; coastsegtype = []
    coastfile = 'gshhs_'+resolution+'.txt'
    countryfile = 'countries_'+resolution+'.txt'
    statefile = 'states_'+resolution+'.txt'
    riverfile = 'rivers_'+resolution+'.txt'
    
    poly, polymeta = get_coast_polygons(coastfile)
    f = open('../lib/matplotlib/toolkits/basemap/data/gshhs_'+resolution+'.dat','wb')
    f2 = open('../lib/matplotlib/toolkits/basemap/data/gshhsmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        typ = pm[0]; area = pm[1]; south = pm[2]; north = pm[3]; npts = pm[4]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (typ, area, npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()
    
    poly, polymeta = get_boundary_lines(countryfile)
    f = open('../lib/matplotlib/toolkits/basemap/data/countries_'+resolution+'.dat','wb')
    f2 = open('../lib/matplotlib/toolkits/basemap/data/countriesmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        south,north,npts = pm[:]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (-1,-1,npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()
    
    poly, polymeta = get_boundary_lines(statefile)
    f = open('../lib/matplotlib/toolkits/basemap/data/states_'+resolution+'.dat','wb')
    f2 = open('../lib/matplotlib/toolkits/basemap/data/statesmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        south,north,npts = pm[:]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (-1,-1,npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()
    
    poly, polymeta = get_boundary_lines(riverfile)
    f = open('../lib/matplotlib/toolkits/basemap/data/rivers_'+resolution+'.dat','wb')
    f2 = open('../lib/matplotlib/toolkits/basemap/data/riversmeta_'+resolution+'.dat','w')
    offset = 0
    for p,pm in zip(poly,polymeta):
        bstring = p.tostring()
        f.write(bstring)
        south,north,npts = pm[:]
        f2.write('%s %s %s %9.5f %9.5f %s %s\n' % (-1,-1,npts, south, north, offset, len(bstring)))
        offset = offset + len(bstring)
    f.close()
    f2.close()
