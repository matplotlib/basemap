# geodesicline.py
#
# This is a rather literal translation of the
# GeographicLib::GeodesicLine class to python.  See the documentation
# for the C++ class for more information at
#
#    http://sourceforge.net/html/annotated.html
#
# The algorithms are derived in
#
#    Charles F. F. Karney,
#    Geodesics on an ellipsoid of revolution, Feb. 2011,
#    http://arxiv.org/abs/1102.1215
#    errata: http://sourceforge.net/geod-errata.html
#
#    Charles F. F. Karney,
#    Algorithms for geodesics, Sept. 2011,
#    http://arxiv.org/abs/1109.4448
#
# Copyright (c) Charles Karney (2011) <charles@karney.com> and licensed
# under the MIT/X11 License.  For more information, see
# http://sourceforge.net/
#
# $Id: 3ce44a02ffcfde1dcd1f660c8fcfa11d024bb4c5 $
######################################################################

import math
from mpl_toolkits.basemap.geomath import Math
from mpl_toolkits.basemap.geodesiccapability import GeodesicCapability

class GeodesicLine(object):
  """Points on a geodesic path"""

  def __init__(self, geod, lat1, lon1, azi1, caps = GeodesicCapability.ALL):
    from pyproj.geodesic import Geodesic
    self._a = geod._a
    self._f = geod._f
    self._b = geod._b
    self._c2 = geod._c2
    self._f1 = geod._f1
    self._caps = caps | Geodesic.LATITUDE | Geodesic.AZIMUTH

    azi1 = Geodesic.AngNormalize(azi1)
    # Guard against underflow in salp0
    azi1 = Geodesic.AngRound(azi1)
    lon1 = Geodesic.AngNormalize(lon1)
    self._lat1 = lat1
    self._lon1 = lon1
    self._azi1 = azi1
    # alp1 is in [0, pi]
    alp1 = azi1 * Math.degree
    # Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
    # problems directly than to skirt them.
    if azi1 == -180:
      self._salp1 = 0
    else:
      self._salp1 = math.sin(alp1)
    if abs(azi1) == 90:
      self._calp1 = 0
    else:
      self._calp1 = math.cos(alp1)
    # real cbet1, sbet1, phi
    phi = lat1 * Math.degree
    # Ensure cbet1 = +epsilon at poles
    sbet1 = self._f1 * math.sin(phi)
    if abs(lat1) == 90:
      cbet1 = Geodesic.tiny_
    else:
      cbet1 = math.cos(phi)
    sbet1, cbet1 = Geodesic.SinCosNorm(sbet1, cbet1)

    # Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
    self._salp0 = self._salp1 * cbet1 # alp0 in [0, pi/2 - |bet1|]
    # Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
    # is slightly better (consider the case salp1 = 0).
    self._calp0 = math.hypot(self._calp1, self._salp1 * sbet1)
    # Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
    # sig = 0 is nearest northward crossing of equator.
    # With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
    # With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
    # With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
    # Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
    # With alp0 in (0, pi/2], quadrants for sig and omg coincide.
    # No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
    # With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
    self._ssig1 = sbet1; self._somg1 = self._salp0 * sbet1
    if sbet1 != 0 or self._calp1 != 0:
      self._csig1 = self._comg1 = cbet1 * self._calp1
    else:
      self._csig1 = self._comg1 = 1.0
    # sig1 in (-pi, pi]
    self._ssig1, self._csig1 = Geodesic.SinCosNorm(self._ssig1, self._csig1)
    self._somg1, self._comg1 = Geodesic.SinCosNorm(self._somg1, self._comg1)

    self._k2 = Math.sq(self._calp0) * geod._ep2
    eps = self._k2 / (2 * (1 + math.sqrt(1 + self._k2)) + self._k2)

    if self._caps & Geodesic.CAP_C1:
      self._A1m1 = Geodesic.A1m1f(eps)
      self._C1a = list(range(Geodesic.nC1_ + 1))
      Geodesic.C1f(eps, self._C1a)
      self._B11 = Geodesic.SinCosSeries(
        True, self._ssig1, self._csig1, self._C1a, Geodesic.nC1_)
      s = math.sin(self._B11); c = math.cos(self._B11)
      # tau1 = sig1 + B11
      self._stau1 = self._ssig1 * c + self._csig1 * s
      self._ctau1 = self._csig1 * c - self._ssig1 * s
      # Not necessary because C1pa reverts C1a
      #    _B11 = -SinCosSeries(true, _stau1, _ctau1, _C1pa, nC1p_)

    if self._caps & Geodesic.CAP_C1p:
      self._C1pa = list(range(Geodesic.nC1p_ + 1))
      Geodesic.C1pf(eps, self._C1pa)

    if self._caps & Geodesic.CAP_C2:
      self._A2m1 =  Geodesic.A2m1f(eps)
      self._C2a = list(range(Geodesic.nC2_ + 1))
      Geodesic.C2f(eps, self._C2a)
      self._B21 = Geodesic.SinCosSeries(
        True, self._ssig1, self._csig1, self._C2a, Geodesic.nC2_)

    if self._caps & Geodesic.CAP_C3:
      self._C3a = list(range(Geodesic.nC3_))
      geod.C3f(eps, self._C3a)
      self._A3c = -self._f * self._salp0 * geod.A3f(eps)
      self._B31 = Geodesic.SinCosSeries(
        True, self._ssig1, self._csig1, self._C3a, Geodesic.nC3_-1)

    if self._caps & Geodesic.CAP_C4:
      self._C4a = list(range(Geodesic.nC4_))
      geod.C4f(self._k2, self._C4a)
      # Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
      self._A4 = Math.sq(self._a) * self._calp0 * self._salp0 * geod._e2
      self._B41 = Geodesic.SinCosSeries(
        False, self._ssig1, self._csig1, self._C4a, Geodesic.nC4_)

  # return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12
  def GenPosition(self, arcmode, s12_a12, outmask):

    from pyproj.geodesic import Geodesic
    a12 = lat2 = lon2 = azi2 = s12 = m12 = M12 = M21 = S12 = Math.nan
    outmask &= self._caps & Geodesic.OUT_ALL
    if not (arcmode or (self._caps & Geodesic.DISTANCE_IN & Geodesic.OUT_ALL)):
      # Uninitialized or impossible distance calculation requested
      return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12

    # Avoid warning about uninitialized B12.
    B12 = 0; AB1 = 0
    if arcmode:
      # Interpret s12_a12 as spherical arc length
      sig12 = s12_a12 * Math.degree
      s12a = abs(s12_a12)
      s12a -= 180 * math.floor(s12a / 180)
      if s12a == 0:
        ssig12 = 0.0
      else:
        ssig12 = math.sin(sig12)
      if s12a == 90:
        csig12 = 0.0
      else:
        csig12 = math.cos(sig12)
    else:
      # Interpret s12_a12 as distance
      tau12 = s12_a12 / (self._b * (1 + self._A1m1))
      s = math.sin(tau12); c = math.cos(tau12)
      # tau2 = tau1 + tau12
      B12 = - Geodesic.SinCosSeries(True, self._stau1 * c + self._ctau1 * s,
                                     self._ctau1 * c - self._stau1 * s,
                                     self._C1pa, Geodesic.nC1p_)
      sig12 = tau12 - (B12 - self._B11)
      ssig12 = math.sin(sig12); csig12 = math.cos(sig12)

    # real omg12, lam12, lon12
    # real ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2
    # sig2 = sig1 + sig12
    ssig2 = self._ssig1 * csig12 + self._csig1 * ssig12
    csig2 = self._csig1 * csig12 - self._ssig1 * ssig12
    if outmask & (
      Geodesic.DISTANCE | Geodesic.REDUCEDLENGTH | Geodesic.GEODESICSCALE):
      if arcmode:
        B12 = Geodesic.SinCosSeries(True, ssig2, csig2,
                                    self._C1a, Geodesic.nC1_)
      AB1 = (1 + self._A1m1) * (B12 - self._B11)
    # sin(bet2) = cos(alp0) * sin(sig2)
    sbet2 = self._calp0 * ssig2
    # Alt: cbet2 = hypot(csig2, salp0 * ssig2)
    cbet2 = math.hypot(self._salp0, self._calp0 * csig2)
    if cbet2 == 0:
      # I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
      cbet2 = csig2 = Geodesic.tiny_
    # tan(omg2) = sin(alp0) * tan(sig2)
    somg2 = self._salp0 * ssig2; comg2 = csig2  # No need to normalize
    # tan(alp0) = cos(sig2)*tan(alp2)
    salp2 = self._salp0; calp2 = self._calp0 * csig2 # No need to normalize
    # omg12 = omg2 - omg1
    omg12 = math.atan2(somg2 * self._comg1 - comg2 * self._somg1,
                  comg2 * self._comg1 + somg2 * self._somg1)

    if outmask & Geodesic.DISTANCE:
      if arcmode:
        s12 = self._b * ((1 + self._A1m1) * sig12 + AB1)
      else:
        s12 = s12_a12

    if outmask & Geodesic.LONGITUDE:
      lam12 = omg12 + self._A3c * (
        sig12 + (Geodesic.SinCosSeries(True, ssig2, csig2,
                                       self._C3a, Geodesic.nC3_-1)
                 - self._B31))
      lon12 = lam12 / Math.degree
      # Can't use AngNormalize because longitude might have wrapped multiple
      # times.
      lon12 = lon12 - 360 * math.floor(lon12/360 + 0.5)
      lon2 = Geodesic.AngNormalize(self._lon1 + lon12)

    if outmask & Geodesic.LATITUDE:
      lat2 = math.atan2(sbet2, self._f1 * cbet2) / Math.degree

    if outmask & Geodesic.AZIMUTH:
      # minus signs give range [-180, 180). 0- converts -0 to +0.
      azi2 = 0 - math.atan2(-salp2, calp2) / Math.degree

    if outmask & (Geodesic.REDUCEDLENGTH | Geodesic.GEODESICSCALE):
      ssig1sq = Math.sq(self._ssig1)
      ssig2sq = Math.sq( ssig2)
      w1 = math.sqrt(1 + self._k2 * ssig1sq)
      w2 = math.sqrt(1 + self._k2 * ssig2sq)
      B22 = Geodesic.SinCosSeries(True, ssig2, csig2, self._C2a, Geodesic.nC2_)
      AB2 = (1 + self._A2m1) * (B22 - self._B21)
      J12 = (self._A1m1 - self._A2m1) * sig12 + (AB1 - AB2)
      if outmask & Geodesic.REDUCEDLENGTH:
        # Add parens around (_csig1 * ssig2) and (_ssig1 * csig2) to ensure
        # accurate cancellation in the case of coincident points.
        m12 = self._b * ((w2 * (self._csig1 * ssig2) -
                          w1 * (self._ssig1 * csig2))
                  - self._csig1 * csig2 * J12)
      if outmask & Geodesic.GEODESICSCALE:
        M12 = csig12 + (self._k2 * (ssig2sq - ssig1sq) *  ssig2 / (w1 + w2)
                        - csig2 * J12) * self._ssig1 / w1
        M21 = csig12 - (self._k2 * (ssig2sq - ssig1sq) * self._ssig1 / (w1 + w2)
                        - self._csig1 * J12) * ssig2 / w2

    if outmask & Geodesic.AREA:
      B42 = Geodesic.SinCosSeries(False, ssig2, csig2, self._C4a, Geodesic.nC4_)
      # real salp12, calp12
      if self._calp0 == 0 or self._salp0 == 0:
        # alp12 = alp2 - alp1, used in atan2 so no need to normalized
        salp12 = salp2 * self._calp1 - calp2 * self._salp1
        calp12 = calp2 * self._calp1 + salp2 * self._salp1
        # The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
        # salp12 = -0 and alp12 = -180.  However this depends on the sign being
        # attached to 0 correctly.  The following ensures the correct behavior.
        if salp12 == 0 and calp12 < 0:
          salp12 = Geodesic.tiny_ * self._calp1
          calp12 = -1
      else:
        # tan(alp) = tan(alp0) * sec(sig)
        # tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
        # = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
        # If csig12 > 0, write
        #   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
        # else
        #   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
        # No need to normalize
        if csig12 > 0:
          salp12 = self._calp0 * self._salp0 * (
            ssig12 * (self._csig1 * ssig12 / (1 + csig12) + self._ssig1))
        else:
          salp12 = self._calp0 * self._salp0 * (
            self._csig1 * (1 - csig12) + ssig12 * self._ssig1)
        calp12 = (Math.sq(self._salp0) +
                  Math.sq(self._calp0) * self._csig1 * csig2)
      S12 = self._c2 * math.atan2(salp12, calp12) + self._A4 * (B42 - self._B41)

    if arcmode:
      a12 = s12_a12
    else:
      a12 = sig12 / Math.degree
    return a12, lat2, lon2, azi2, s12, m12, M12, M21, S12

  def Position(self, s12,
               outmask = GeodesicCapability.LATITUDE |
               GeodesicCapability.LONGITUDE | GeodesicCapability.AZIMUTH):
    """
    Return the point a distance s12 along the geodesic line.  Return
    a dictionary with (some) of the following entries:

      lat1 latitude of point 1
      lon1 longitude of point 1
      azi1 azimuth of line at point 1
      lat2 latitude of point 2
      lon2 longitude of point 2
      azi2 azimuth of line at point 2
      s12 distance from 1 to 2
      a12 arc length on auxiliary sphere from 1 to 2
      m12 reduced length of geodesic
      M12 geodesic scale 2 relative to 1
      M21 geodesic scale 1 relative to 2
      S12 area between geodesic and equator

    outmask determines which fields get included and if outmask is
    omitted, then only the basic geodesic fields are computed.  The mask
    is an or'ed combination of the following values

      Geodesic.LATITUDE
      Geodesic.LONGITUDE
      Geodesic.AZIMUTH
      Geodesic.DISTANCE
      Geodesic.REDUCEDLENGTH
      Geodesic.GEODESICSCALE
      Geodesic.AREA
      Geodesic.ALL
    """

    from pyproj.geodesic import Geodesic
    Geodesic.CheckDistance(s12)
    result = {'lat1': self._lat1, 'lon1': self._lon1, 'azi1': self._azi1,
              's12': s12}
    a12, lat2, lon2, azi2, s12, m12, M12, M21, S12 = self.GenPosition(
      False, s12, outmask)
    outmask &= Geodesic.OUT_ALL
    result['a12'] = a12
    if outmask & Geodesic.LATITUDE: result['lat2'] = lat2
    if outmask & Geodesic.LONGITUDE: result['lon2'] = lon2
    if outmask & Geodesic.AZIMUTH: result['azi2'] = azi2
    if outmask & Geodesic.REDUCEDLENGTH: result['m12'] = m12
    if outmask & Geodesic.GEODESICSCALE:
      result['M12'] = M12; result['M21'] = M21
    if outmask & Geodesic.AREA: result['S12'] = S12
    return result
