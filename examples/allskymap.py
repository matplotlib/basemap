from __future__ import unicode_literals
"""
AllSkyMap is a subclass of Basemap, specialized for handling common plotting
tasks for celestial data.

It is essentially equivalent to using Basemap with full-sphere projections
(e.g., 'hammer' or 'moll') and the `celestial` keyword set to `True`, but
it adds a few new methods:

* label_meridians for, well, labeling meridians with their longitude values;

* geodesic, a replacement for Basemap.drawgreatcircle, that can correctly
  handle geodesics that cross the limb of the map, and providing the user
  easy control over clipping (which affects thick lines at or near the limb);
  
* tissot, which overrides Basemap.tissot, correctly handling geodesics that
  cross the limb of the map.

Created Jan 2011 by Tom Loredo, based on Jeff Whitaker's code in Basemap's
__init__.py module.
"""

from numpy import *
import matplotlib.pyplot as pl
from matplotlib.pyplot import *
from mpl_toolkits.basemap import Basemap, pyproj
from mpl_toolkits.basemap.pyproj import Geod

__all__ = ['AllSkyMap']

def angle_symbol(angle, round_to=1.0):
    """
    Return a string representing an angle, rounded and with a degree symbol.
    
    This is adapted from code in mpl's projections.geo module.
    """
    value = np.round(angle / round_to) * round_to
    if pl.rcParams['text.usetex'] and not pl.rcParams['text.latex.unicode']:
        return r'$%0.0f^\circ$' % value
    else:
        return '%0.0f\N{DEGREE SIGN}' % value


class AllSkyMap(Basemap):
    """
    AllSkyMap is a subclass of Basemap, specialized for handling common plotting
    tasks for celestial data.
    
    It is essentially equivalent to using Basemap with full-sphere projections
    (e.g., 'hammer' or 'moll') and the `celestial` keyword set to `True`, but
    it adds a few new methods:
    
    * label_meridians for, well, labeling meridians with their longitude values;
    
    * geodesic, a replacement for Basemap.drawgreatcircle, that can correctly
      handle geodesics that cross the limb of the map, and providing the user
      easy control over clipping (which affects thick lines at or near the
      limb);
      
    * tissot, which overrides Basemap.tissot, correctly handling geodesics that
      cross the limb of the map.
    """

    # Longitudes corresponding to east and west edges, reflecting the
    # convention that 180 deg is the eastern edge, according to basemap's 
    # underlying projections:
    east_lon = 180.
    west_lon = 180.+1.e-10

    def __init__(self, 
                       projection='hammer',
                       lat_0=0., lon_0=0.,
                       suppress_ticks=True,
                       boundinglat=None,
                       fix_aspect=True,
                       anchor=str('C'),
                       ax=None):

        if projection != 'hammer' and projection !='moll':
            raise ValueError('Only hammer and moll projections supported!')

        # Use Basemap's init, enforcing the values of many parameters that
        # aren't used or whose Basemap defaults would not be altered for all-sky
        # celestial maps.
        Basemap.__init__(self, llcrnrlon=None, llcrnrlat=None,
                       urcrnrlon=None, urcrnrlat=None,
                       llcrnrx=None, llcrnry=None,
                       urcrnrx=None, urcrnry=None,
                       width=None, height=None,
                       projection=projection, resolution=None,
                       area_thresh=None, rsphere=1.,
                       lat_ts=None,
                       lat_1=None, lat_2=None,
                       lat_0=lat_0, lon_0=lon_0,
                       suppress_ticks=suppress_ticks,
                       satellite_height=1.,
                       boundinglat=None,
                       fix_aspect=True,
                       anchor=anchor,
                       celestial=True,
                       ax=ax)

        # Keep a local ref to lon_0 for hemisphere checking.
        self._lon_0 = self.projparams['lon_0']
        self._limb = None

    def drawmapboundary(self,color='k',linewidth=1.0,fill_color=None,\
                        zorder=None,ax=None):
        """
        draw boundary around map projection region, optionally
        filling interior of region.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        line width for boundary (default 1.)
        color            color of boundary line (default black)
        fill_color       fill the map region background with this
                         color (default is no fill or fill with axis
                         background color).
        zorder           sets the zorder for filling map background
                         (default 0).
        ax               axes instance to use
                         (default None, use default axes instance).
        ==============   ====================================================

        returns matplotlib.collections.PatchCollection representing map boundary.
        """
        # Just call the base class version, but keep a copy of the limb
        # polygon for clipping.
        self._limb = Basemap.drawmapboundary(self, color=color,
            linewidth=linewidth, fill_color=fill_color, zorder=zorder, ax=ax)
        return self._limb

    def label_meridians(self, lons, fontsize=10, valign='bottom', vnudge=0,
                        halign='center', hnudge=0):
        """
        Label meridians with their longitude values in degrees.
        
        This labels meridians with negative longitude l with the value 360-l;
        for maps in celestial orientation, this means meridians to the right
        of the central meridian are labeled from 360 to 180 (left to right).
        
        `vnudge` and `hnudge` specify amounts in degress to nudge the labels
        from their default placements, vertically and horizontally.  This
        values obey the map orientation, so to nudge to the right, use a
        negative `hnudge` value.
        """
        # Run through (lon, lat) pairs, with lat=0 in each pair.
        lats = len(lons)*[0.]
        for lon,lat in zip(lons, lats):
            x, y = self(lon+hnudge, lat+vnudge)
            if lon < 0:
                lon_lbl = 360 + lon
            else:
                lon_lbl = lon
            pl.text(x, y, angle_symbol(lon_lbl), fontsize=fontsize,
                    verticalalignment=valign,
                    horizontalalignment=halign)

    def east_hem(self, lon):
        """
        Return True if lon is in the eastern hemisphere of the map wrt lon_0.
        """
        if (lon-self._lon_0) % 360. <= self.east_lon:
            return True
        else:
            return False

    def geodesic(self, lon1, lat1, lon2, lat2, del_s=.01, clip=True, **kwargs):
        """
        Plot a geodesic curve from (lon1, lat1) to (lon2, lat2), with
        points separated by arc length del_s.  Return a list of Line2D
        instances for the curves comprising the geodesic.  If the geodesic does
        not cross the map limb, there will be only a single curve; if it
        crosses the limb, there will be two curves.
        """
        
        # TODO:  Perhaps return a single Line2D instance when there is only a
        # single segment, and a list of segments only when there are two segs?

        # TODO:  Check the units of del_s.
        
        # This is based on Basemap.drawgreatcircle (which draws an *arc* of a
        # great circle), but addresses a limitation of that method, supporting
        # geodesics that cross the map boundary by breaking them into two
        # segments, one in the eastern hemisphere and the other in the western.
        gc = pyproj.Geod(a=self.rmajor,b=self.rminor)
        az12,az21,dist = gc.inv(lon1,lat1,lon2,lat2)
        npoints = int((dist+0.5**del_s)/del_s)
        # Calculate lon & lat for points on the arc.
        lonlats = gc.npts(lon1,lat1,lon2,lat2,npoints)
        lons = [lon1]; lats = [lat1]
        for lon, lat in lonlats:
            lons.append(lon)
            lats.append(lat)
        lons.append(lon2); lats.append(lat2)
        # Break the arc into segments as needed, when there is a longitudinal
        # hemisphere crossing.
        segs = []
        seg_lons, seg_lats = [lon1], [lat1]
        cur_hem = self.east_hem(lon1)
        for lon, lat in zip(lons[1:], lats[1:]):
            if self.east_hem(lon) == cur_hem:
                seg_lons.append(lon)
                seg_lats.append(lat)
            else:
                # We should interpolate a new pt at the boundary, but in
                # the mean time just rely on the step size being small.
                segs.append( (seg_lons, seg_lats) )
                seg_lons, seg_lats = [lon], [lat]
                cur_hem = not cur_hem
        segs.append( (seg_lons, seg_lats) )
        # Plot each segment; return a list of the mpl lines.
        lines = []
        for lons, lats in segs:
            x, y = self(lons, lats)
            if clip and self._limb:
                line = plot(x, y, clip_path=self._limb, **kwargs)[0]
            else:
                line = plot(x, y, **kwargs)[0]
            lines.append(line)
        # If there are multiple segments and no color args, reconcile the
        # colors, which mpl will have autoset to different values.
        # *** Does this screw up mpl's color set sequence for later lines?
        if 'c' not in kwargs or 'color' in kwargs:
            if len(lines) > 1:
                c1 = lines[0].get_color()
                for line in lines[1:]:
                    line.set_color(c1)
        return lines

    def tissot(self,lon_0,lat_0,radius_deg,npts,ax=None,**kwargs):
        """
        Draw a polygon centered at ``lon_0,lat_0``.  The polygon
        approximates a circle on the surface of the earth with radius
        ``radius_deg`` degrees latitude along longitude ``lon_0``,
        made up of ``npts`` vertices.
        
        The polygon represents a Tissot's indicatrix
        (http://en.wikipedia.org/wiki/Tissot's_Indicatrix),
        which when drawn on a map shows the distortion inherent in the map
        projection.  Tissots can be used to display azimuthally symmetric
        directional uncertainties ("error circles").

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.patches.Polygon.

        returns a list of matplotlib.patches.Polygon objects, with two polygons
        when the tissot crosses the limb, and just one polygon otherwise.
        """
        
        # TODO:  Just return the polygon (not a list) when there is only one
        # polygon?  Or stick with the list for consistency?
        
        # This is based on Basemap.tissot, but addresses a limitation of that
        # method by handling tissots that cross the limb of the map by finding
        # separate polygons in the eastern and western hemispheres comprising
        # the tissot.
        ax = kwargs.pop('ax', None) or self._check_ax()
        g = pyproj.Geod(a=self.rmajor,b=self.rminor)
        az12,az21,dist = g.inv(lon_0,lat_0,lon_0,lat_0+radius_deg)
        start_hem = self.east_hem(lon_0)
        segs1 = [self(lon_0,lat_0+radius_deg)]
        over, segs2 = [], []
        delaz = 360./npts
        az = az12
        last_lon = lon_0
        # Note adjacent and opposite edge longitudes, in case the tissot
        # runs over the edge.
        if start_hem:  # eastern case
            adj_lon = self.east_lon
            opp_lon = self.west_lon
        else:
            adj_lon = self.west_lon
            opp_lon = self.east_lon
        for n in range(npts):
            az = az+delaz
            # skip segments along equator (Geod can't handle equatorial arcs)
            if np.allclose(0.,lat_0) and (np.allclose(90.,az) or np.allclose(270.,az)):
                continue
            else:
                lon, lat, az21 = g.fwd(lon_0, lat_0, az, dist)
            # If in the starting hemisphere, add to 1st polygon seg list.
            if self.east_hem(lon) == start_hem:
                x, y = self(lon, lat)
                # Add segment if it is in the map projection region.
                if x < 1.e20 and y < 1.e20:
                    segs1.append( (x, y) )
                    last_lon = lon
            # Otherwise, we cross hemispheres.
            else:
                # Trace the edge of each hemisphere.
                x, y = self(adj_lon, lat)
                if x < 1.e20 and y < 1.e20:
                    segs1.append( (x, y) )
                    # We presume if adj projection is okay, opposite is.
                    segs2.append( self(opp_lon, lat) )
                # Also store the overlap in the opposite hemisphere.
                x, y = self(lon, lat)
                if x < 1.e20 and y < 1.e20:
                    over.append( (x, y) )
                    last_lon = lon
        poly1 = Polygon(segs1, **kwargs)
        ax.add_patch(poly1)
        if segs2:
            over.reverse()
            segs2.extend(over)
            poly2 = Polygon(segs2, **kwargs)
            ax.add_patch(poly2)
            return [poly1, poly2]
        else:
            return [poly1]


if __name__ == '__main__':

    # Note that Hammer & Mollweide projections enforce a 2:1 aspect ratio.
    # Use figure size good for a 2:1 plot.
    fig = figure(figsize=(12,6))
    
    # Set up the projection and draw a grid.
    map = AllSkyMap(projection='hammer')
    # Save the bounding limb to use as a clip path later.
    limb = map.drawmapboundary(fill_color='white')
    map.drawparallels(np.arange(-75,76,15), linewidth=0.5, dashes=[1,2],
        labels=[1,0,0,0], fontsize=9)
    map.drawmeridians(np.arange(-150,151,30), linewidth=0.5, dashes=[1,2])
    
    # Label a subset of meridians.
    lons = np.arange(-150,151,30)
    map.label_meridians(lons, fontsize=9, vnudge=1,
                    halign='left', hnudge=-1)  # hnudge<0 shifts to right
    
    # x, y limits are [0, 4*rt2], [0, 2*rt2].
    rt2 = sqrt(2)

    # Draw a slanted green line crossing the map limb.
    line = plot([rt2,0], [rt2,2*rt2], 'g-')

    # Draw a slanted magenta line crossing the map limb but clipped.
    line = plot([rt2+.1,0+.1], [rt2,2*rt2], 'm-', clip_path=limb)
    
    # Draw some geodesics.
    # First a transparent thick blue geodesic crossing the limb but not clipped,
    # overlayed by a thinner red geodesic that is clipped (by default), to
    # illustrate the effect of clipping.
    lines = map.geodesic(120, 30, 240, 60, clip=False, c='b', lw=7, alpha=.5)
    lines = map.geodesic(240, 60, 120, 30, c='r', lw=3, alpha=.5)

    # Next two large limb-crossing geodesics with the same path, but rendered
    # in opposite directions, one transparent blue, the other transparent
    # yellow.  They should be right on top of each other, giving a greenish
    # brown hue.
    lines = map.geodesic(240, -60, 120, 30, c='b', lw=2, alpha=.5)
    lines = map.geodesic(120, 30, 240, -60, c='y', lw=2, alpha=.5)

    # What happens if a geodesic is given coordinates spanning more than
    # a single rotation?  Not sure what to expect, but it shoots off the
    # map (clipped here).  Perhaps we should ensure lons are in [0, 360].
    #lines = map.geodesic(120, 20, 240+360, 50, del_s=.2, c='g')
    
    # Two tissots fully within the limb.
    poly = map.tissot(60, -15, 10, 100)
    poly = map.tissot(280, 60, 10, 100)
    #poly = map.tissot(90, -85, 10, 100)
    
    # Limb-spanning tissots in each quadrant.
    # lower left:
    poly = map.tissot(170, -60, 15, 100)
    # upper left:
    poly = map.tissot(175, 70, 15, 100)
    # upper right (note negative longitude):
    poly = map.tissot(-175, 30, 15, 100, color='r', alpha=.6)
    # lower right:
    poly = map.tissot(185, -40, 10, 100)

    # Plot the tissot centers as "+" symbols.  Note the top left symbol
    # would cross the limb without the clip_path argument; this might be
    # desired to enhance visibility.
    lons = [170, 175, -175, 185]
    lats = [-60, 70, 30, -40]
    x, y = map(lons, lats)
    map.scatter(x, y, s=40, marker='+', linewidths=1, edgecolors='g',
        facecolors='none', clip_path=limb, zorder=10)  # hi zorder -> top
    
    title('AllSkyMap demo:  Clipped lines, markers, geodesics, tissots')
    show()
