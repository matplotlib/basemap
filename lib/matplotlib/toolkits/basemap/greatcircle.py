import math
import numarray as N

class GreatCircle:
    """
    from Ed Williams' 'Aviation Formulary'
    (http://williams.best.vwh.net/avform.htm)
    """
    def __init__(self,lon1,lat1,lon2,lat2):
        """
        Define a great circle by specifying:
        lon1 - starting longitude of great circle
        lat1 - starting latitude
        lon2 - ending longitude
        lat2 - ending latitude
        All must be given in degrees.

        Instance variables:
        distance - distance along great circle in radians.
        lon1,lat1,lon2,lat2 - start and end points (in radians).
        antipodal - True if start and end points are antipodal.
        """
        # convert to radians from degrees.
        lat1 = math.radians(lat1)
        lon1 = math.radians(lon1)
        lat2 = math.radians(lat2)
        lon2 = math.radians(lon2)
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        # distance along great circle in radians.
        self.distance = 2.*math.asin(math.sqrt((math.sin((lat1-lat2)/2))**2+\
        math.cos(lat1)*math.cos(lat2)*(math.sin((lon1-lon2)/2))**2))
        # check to see if points are antipodal (if so, route is undefined).
        if self.distance == math.pi:
            self.antipodal = True
        else:
            self.antipodal = False

    def points(self,npoints):
        """
        compute arrays of npoints equally spaced
        intermediate points along the great circle.

        input parameter npoints is the number of points
        to compute.

        Returns lons, lats (lists with longitudes and latitudes
        of intermediate points in degrees).

        For example npoints=10 will return arrays lons,lats of 10
        equally spaced points along the great circle.
        """
        # can't do it if endpoints of great circle are antipodal, since
        # route is undefined.
        if self.antipodal:
            raise ValueError,'cannot compute intermediate points on a great circle whose endpoints are antipodal'
        d = self.distance
        delta = 1.0/(npoints-1)
        f = delta*N.arange(npoints)
        lat1 = self.lat1
        lat2 = self.lat2
        lon1 = self.lon1
        lon2 = self.lon2
        A = N.sin((1-f)*d)/math.sin(d)
        B = N.sin(f*d)/math.sin(d)
        x = A*math.cos(lat1)*math.cos(lon1)+B*math.cos(lat2)*math.cos(lon2)
        y = A*math.cos(lat1)*math.sin(lon1)+B*math.cos(lat2)*math.sin(lon2)
        z = A*math.sin(lat1)               +B*math.sin(lat2)
        lats=N.arctan2(z,N.sqrt(x**2+y**2))
        lons=N.arctan2(y,x)
        lons = map(math.degrees,lons.tolist())
        lats = map(math.degrees,lats.tolist())
        return lons,lats

if __name__ == '__main__':
    lon1 = float(raw_input('enter longitude of start point:'))
    lat1 = float(raw_input('enter latitude of start point:'))
    lon2 = float(raw_input('enter longitude of end point:'))
    lat2 = float(raw_input('enter latitude of end point:'))
    gc = GreatCircle(lon1,lat1,lon2,lat2)
    print 'distance along great circle (in radians) = ',gc.distance
    print 'lon/lat for 10 equally spaced points along great circle:'
    lons,lats = gc.points(10)
    for lon,lat in zip(lons,lats):
        print lon,lat
