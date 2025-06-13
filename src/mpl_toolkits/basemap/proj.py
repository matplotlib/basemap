"""Module with :mod:`Proj` class for cartographic transformations."""
from __future__ import division

import math
import numpy as np
import pyproj


_dg2rad = math.radians(1.)
_rad2dg = math.degrees(1.)

_cylproj = ["cyl", "merc", "mill", "gall"]
_pseudocyl = ["moll", "kav7", "eck4", "robin", "sinu", "mbtfpq", "vandg", "hammer"]


class Proj(object):  # pylint: disable=too-many-instance-attributes
    """Perform cartographic transformations using :mod:`pyproj`.

    The cartographic transformations (from longitude/latitude to native
    map projection coorinates x/y coordinates and vice versa) is done
    with :mod:`pyproj`, a Cython interface to PROJ (https://proj.org/).

    * __init__ method sets up projection information.
    * __call__ method compute transformations.

    See docstrings for __init__ and __call__ for details."""

    def __init__(self, projparams,  # pylint: disable=too-many-arguments
                 llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, urcrnrislatlon=True):
        """Initialise a :mod:`Proj` instance.

        Input `projparams` is a dictionary with PROJ map projection
        control parameters given as key/value pairs. See the PROJ
        documentation (https://proj.org/) for details.

        llcrnrlon,llcrnrlat are lon and lat (in degrees) of lower
        left hand corner of projection region.

        urcrnrlon,urcrnrlat are lon and lat (in degrees) of upper
        right hand corner of projection region if urcrnrislatlon=True
        (default). Otherwise, urcrnrlon,urcrnrlat are x,y in projection
        coordinates (units meters), assuming the lower left corner is
        x=0,y=0.
        """

        # pylint: disable=too-many-statements,too-many-branches
        self.projparams = projparams
        self.projection = projparams["proj"]

        # rmajor is the semi-major axis.
        # rminor is the semi-minor axis.
        # esq is eccentricity squared.
        try:
            self.rmajor = projparams["a"]
            self.rminor = projparams["b"]
        except:  # noqa: E722  # pylint: disable=bare-except
            try:
                self.rmajor = projparams["R"]
            except:  # noqa: E722  # pylint: disable=bare-except
                self.rmajor = projparams["bR_a"]
            self.rminor = self.rmajor
        if self.rmajor == self.rminor:
            self.ellipsoid = False
        else:
            self.ellipsoid = True
        self.flattening = (self.rmajor - self.rminor) / self.rmajor
        self.esq = (self.rmajor**2 - self.rminor**2) / self.rmajor**2

        self.llcrnrlon = llcrnrlon
        self.llcrnrlat = llcrnrlat
        if self.projection == "cyl":
            llcrnrx = llcrnrlon
            llcrnry = llcrnrlat
        elif self.projection == "ob_tran":
            self._proj4 = pyproj.Proj(projparams)
            llcrnrx, llcrnry = self(llcrnrlon, llcrnrlat)
            llcrnrx = _rad2dg * llcrnrx
            llcrnry = _rad2dg * llcrnry
            if llcrnrx < 0:
                llcrnrx = llcrnrx + 360
        elif self.projection in "ortho":
            if llcrnrlon == -180 and llcrnrlat == -90 and urcrnrlon == +180 and urcrnrlat == +90:
                self._fulldisk = True
                self._proj4 = pyproj.Proj(projparams)
                llcrnrx = -self.rmajor
                llcrnry = -self.rmajor
                self._width = 0.5 * (self.rmajor + self.rminor)
                self._height = 0.5 * (self.rmajor + self.rminor)
                urcrnrx = -llcrnrx
                urcrnry = -llcrnry
            else:
                self._fulldisk = False
                self._proj4 = pyproj.Proj(projparams)
                llcrnrx, llcrnry = self(llcrnrlon, llcrnrlat)
                if llcrnrx > 1E20 or llcrnry > 1E20:
                    raise ValueError("the lower left corner of the plot "
                                     "is not in the map projection region")
        elif self.projection == "aeqd" and (llcrnrlon == -180 and llcrnrlat == -90 and
                                            urcrnrlon == +180 and urcrnrlat == +90):
            self._fulldisk = True
            self._proj4 = pyproj.Proj(projparams)
            # Raise an exception for ellipsoids - there appears to be a bug
            # in PROJ4 that causes the inverse transform to fail for points
            # more than 90 degrees of arc away from center point for ellipsoids
            # (works fine for spheres) - below is an example
            # from pyproj import Proj
            # p1 = Proj(proj="aeqd", a=6378137.0, b=6356752.3142, lat_0=0, lon_0=0)
            # x, y = p1(91, 0)
            # lon, lat = p1(x, y, inverse=True)  # lon is 89 instead of 91
            if self.ellipsoid:
                raise ValueError(
                    "full disk (whole world) Azimuthal Equidistant projection "
                    "can only be drawn for a perfect sphere")
            llcrnrx = -np.pi * self.rmajor
            llcrnry = -np.pi * self.rmajor
            self._width = -llcrnrx
            self._height = -llcrnry
            urcrnrx = -llcrnrx
            urcrnry = -llcrnry
        elif self.projection == "geos":
            self._proj4 = pyproj.Proj(projparams)
            # Find major and minor axes of ellipse defining map proj region.
            # h is measured from surface of earth at equator.
            h = projparams["h"] + self.rmajor
            # Latitude of horizon on central meridian.
            lonmax = 90. - (180. / np.pi) * np.arcsin(self.rmajor / h)
            # Longitude of horizon on equator.
            latmax = 90. - (180. / np.pi) * np.arcsin(self.rminor / h)
            # Truncate to nearest hundredth of a degree (to make sure
            # they aren't slightly over the horizon).
            latmax = int(100 * latmax) / 100.
            lonmax = int(100 * lonmax) / 100.
            # Width and height of visible projection.
            projtmp = pyproj.Proj(proj="geos", a=self.rmajor, b=self.rminor,
                                  lat_0=0, lon_0=0, h=projparams["h"])
            _x1, y1 = projtmp(0., latmax)  # pylint: disable=unpacking-non-sequence
            x2, _y2 = projtmp(lonmax, 0.)  # pylint: disable=unpacking-non-sequence
            width, height = x2, y1
            self._height = height
            self._width = width
            if llcrnrlon == -180 and llcrnrlat == -90 and urcrnrlon == +180 and urcrnrlat == +90:
                self._fulldisk = True
                llcrnrx = np.negative(width)
                llcrnry = np.negative(height)
                urcrnrx = np.negative(llcrnrx)
                urcrnry = np.negative(llcrnry)
            else:
                self._fulldisk = False
                llcrnrx, llcrnry = self(llcrnrlon, llcrnrlat)
                if llcrnrx > 1E20 or llcrnry > 1E20:
                    raise ValueError("the lower left corner of the plot "
                                     "is not in the map projection region")
        elif self.projection == "nsper":
            self._proj4 = pyproj.Proj(projparams)
            # Find major and minor axes of ellipse defining map proj region.
            # h is measured from surface of earth at equator.
            h = projparams["h"] + self.rmajor
            # Latitude of horizon on central meridian.
            lonmax = 90. - (180. / np.pi) * np.arcsin(self.rmajor / h)
            # Longitude of horizon on equator.
            latmax = 90. - (180. / np.pi) * np.arcsin(self.rmajor / h)
            # Rruncate to nearest hundredth of a degree (to make sure
            # they aren't slightly over the horizon).
            latmax = int(100 * latmax) / 100.
            lonmax = int(100 * lonmax) / 100.
            # Width and height of visible projection.
            projtmp = pyproj.Proj(proj="nsper", a=self.rmajor, b=self.rminor,
                                  lat_0=0, lon_0=0, h=projparams["h"])
            _x1, y1 = projtmp(0., latmax)  # pylint: disable=unpacking-non-sequence
            x2, _y2 = projtmp(lonmax, 0.)  # pylint: disable=unpacking-non-sequence
            width, height = x2, y1
            self._height = height
            self._width = width
            if llcrnrlon == -180 and llcrnrlat == -90 and urcrnrlon == +180 and urcrnrlat == +90:
                self._fulldisk = True
                llcrnrx = np.negative(width)
                llcrnry = np.negative(height)
                urcrnrx = np.negative(llcrnrx)
                urcrnry = np.negative(llcrnry)
            else:
                self._fulldisk = False
                llcrnrx, llcrnry = self(llcrnrlon, llcrnrlat)
                if llcrnrx > 1E20 or llcrnry > 1E20:
                    raise ValueError("the lower left corner of the plot "
                                     "is not in the map projection region")
        elif self.projection in _pseudocyl:
            self._proj4 = pyproj.Proj(projparams)
            _, urcrnry = self(projparams["lon_0"], 90.)
            urcrnrx, _ = self(projparams["lon_0"] + 180. - 1E-10, 0)
            llcrnrx = -urcrnrx
            llcrnry = -urcrnry
            if self.ellipsoid and self.projection in ["kav7", "eck4", "mbtfpq"]:
                msg = "this projection can only be drawn for a perfect sphere"
                raise ValueError(msg)
        else:
            self._proj4 = pyproj.Proj(projparams)
            llcrnrx, llcrnry = self(llcrnrlon, llcrnrlat)
            if self.projection == "aeqd":
                self._fulldisk = False
        # Compute x_0, y_0 so ll corner of domain is x=0,y=0.
        # Note that for "cyl" we have x,y == lon,lat.
        if self.projection != "ob_tran":
            self.projparams["x_0"] = np.negative(llcrnrx)
            self.projparams["y_0"] = np.negative(llcrnry)
        # Reset with x_0, y_0.
        if self.projection not in ["cyl", "ob_tran"]:
            self._proj4 = pyproj.Proj(projparams)
            llcrnry = 0.
            llcrnrx = 0.
        elif self.projection != "ob_tran":
            llcrnrx = llcrnrlon
            llcrnry = llcrnrlat
        if urcrnrislatlon:
            self.urcrnrlon = urcrnrlon
            self.urcrnrlat = urcrnrlat
            if self.projection not in ["ortho", "geos", "nsper", "aeqd"] + _pseudocyl:
                urcrnrx, urcrnry = self(urcrnrlon, urcrnrlat)
                if self.projection == "ob_tran":
                    urcrnrx = _rad2dg * urcrnrx
                    urcrnry = _rad2dg * urcrnry
                    if urcrnrx < 0:
                        urcrnrx = urcrnrx + 360
            elif self.projection in ["ortho", "geos", "nsper", "aeqd"]:
                if self._fulldisk:
                    urcrnrx = 2. * self._width
                    urcrnry = 2. * self._height
                else:
                    urcrnrx, urcrnry = self(urcrnrlon, urcrnrlat)
                    if urcrnrx > 1E20 or urcrnry > 1E20:
                        raise ValueError("the upper right corner of the plot "
                                         "is not in the map projection region")
            elif self.projection in _pseudocyl:
                _, urcrnry = self(projparams["lon_0"], 90.)
                urcrnrx, _ = self(projparams["lon_0"] + 180. - 1E-10, 0)
        else:
            urcrnrx = urcrnrlon
            urcrnry = urcrnrlat
            urcrnrlon, urcrnrlat = self(urcrnrx, urcrnry, inverse=True)
            self.urcrnrlon = urcrnrlon
            self.urcrnrlat = urcrnrlat

        # corners of domain.
        self.llcrnrx = llcrnrx
        self.llcrnry = llcrnry
        self.urcrnrx = urcrnrx
        self.urcrnry = urcrnry
        if urcrnrx > llcrnrx:
            self.xmin = llcrnrx
            self.xmax = urcrnrx
        else:
            self.xmax = llcrnrx
            self.xmin = urcrnrx
        if urcrnry > llcrnry:
            self.ymin = llcrnry
            self.ymax = urcrnry
        else:
            self.ymax = llcrnry
            self.ymin = urcrnry

    def __call__(self, *args, **kw):
        """Perform a cartographic transformation.

        Calling a :mod:`Proj` instance with the arguments lon,lat will
        convert lon/lat (in degrees) to x/y native map projection
        coordinates (in meters). If optional keyword `inverse` is
        True (default is False), the inverse transformation from x/y
        to lon/lat is performed.

        For cylindrical equidistant projection ('cyl'), this
        does nothing (i.e. x,y == lon,lat).

        lon,lat can be either scalar floats or N arrays.
        """

        if len(args) == 1:
            xy = args[0]
            onearray = True
            # Fast return for "cyl" since x,y == lon,lat.
            if self.projection == "cyl":
                return xy
        else:
            x, y = args
            onearray = False
            # Fast return for "cyl" since x,y == lon,lat.
            if self.projection == "cyl":
                return x, y

        inverse = kw.get("inverse", False)
        if onearray:
            outx, outy = self._proj4(*xy.T, inverse=inverse)
            outxy = np.asarray([outx, outy]).T
        else:
            outx, outy = self._proj4(x, y, inverse=inverse)
        if inverse:
            if self.projection in ["merc", "mill", "gall"]:
                if self.projection == "merc":
                    coslat = math.cos(math.radians(self.projparams["lat_ts"]))
                    sinlat = math.sin(math.radians(self.projparams["lat_ts"]))
                else:
                    coslat = 1.
                    sinlat = 0.
                # Radius of curvature of the ellipse perpendicular to
                # the plane of the meridian.
                rcurv = self.rmajor * coslat / math.sqrt(1. - self.esq * sinlat**2)
                if onearray:
                    # pylint: disable=unsupported-assignment-operation
                    outxy[:, 0] = _rad2dg * (xy[:, 0] / rcurv) + self.llcrnrlon
                else:
                    try:
                        # x is a scalar or an array.
                        outx = _rad2dg * (x / rcurv) + self.llcrnrlon
                    except:  # noqa: E722  # pylint: disable=bare-except
                        # x is a sequence.
                        outx = [_rad2dg * (xi / rcurv) + self.llcrnrlon for xi in x]
        else:
            if self.projection in ["merc", "mill", "gall"]:
                if self.projection == "merc":
                    coslat = math.cos(math.radians(self.projparams["lat_ts"]))
                    sinlat = math.sin(math.radians(self.projparams["lat_ts"]))
                else:
                    coslat = 1.
                    sinlat = 0.
                # Radius of curvature of the ellipse perpendicular to
                # the plane of the meridian.
                rcurv = self.rmajor * coslat / math.sqrt(1. - self.esq * sinlat**2)
                if onearray:
                    # pylint: disable=unsupported-assignment-operation
                    outxy[:, 0] = rcurv * _dg2rad * (xy[:, 0] - self.llcrnrlon)
                else:
                    try:
                        # x is a scalar or an array.
                        outx = rcurv * _dg2rad * (x - self.llcrnrlon)
                    except:  # noqa: E722  # pylint: disable=bare-except
                        # x is a sequence.
                        outx = [rcurv * _dg2rad * (xi - self.llcrnrlon) for xi in x]

        if onearray:
            return outxy
        return outx, outy

    def makegrid(self, nx, ny, returnxy=False):
        """Return regular grid in native projection.

        It returns arrays of shape (ny, nx) containing lon,lat
        coordinates of an equally spaced native projection grid.

        If returnxy=True, the x,y values of the grid are returned also.
        """

        dx = (self.urcrnrx - self.llcrnrx) / (nx - 1)
        dy = (self.urcrnry - self.llcrnry) / (ny - 1)
        x = self.llcrnrx + dx * np.indices((ny, nx), np.float32)[1, :, :]
        y = self.llcrnry + dy * np.indices((ny, nx), np.float32)[0, :, :]
        lons, lats = self(x, y, inverse=True)

        if returnxy:
            return lons, lats, x, y
        return lons, lats

    def makegrid3d(self, nx, ny, returnxy=False):
        """Return regular grid in native projection.

        It returns arrays of shape (ny, nx, 2) containing lon,lat
        coordinates of an equally spaced native projection grid.

        if returnxy=True, the x,y values of the grid are returned also.
        """

        dx = (self.urcrnrx - self.llcrnrx) / (nx - 1)
        dy = (self.urcrnry - self.llcrnry) / (ny - 1)
        xy = np.empty((ny, nx, 2), np.float64)

        xy[..., 0] = self.llcrnrx + dx * np.indices((ny, nx), np.float32)[1, :, :]
        xy[..., 1] = self.llcrnry + dy * np.indices((ny, nx), np.float32)[0, :, :]
        lonlat = self(xy, inverse=True)

        if returnxy:
            return lonlat, xy
        return lonlat
