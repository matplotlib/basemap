/*
 * This is a C implementation of the geodesic algorithms described in
 *
 *   C. F. F. Karney
 *   Algorithms for geodesics
 *   J. Geodesy (2012)
 *   http://dx.doi.org/10.1007/s00190-012-0578-z
 *   Addenda: http://geographiclib.sf.net/geod-addenda.html
 *
 * See the comments in geodesic.h for documentation.
 *
 * Copyright (c) Charles Karney (2012) <charles@karney.com> and licensed
 * under the MIT/X11 License.  For more information, see
 * http://geographiclib.sourceforge.net/
 *
 * This file was distributed with GeographicLib 1.27.
 */

#include "geodesic.h"
#include <math.h>

#define GEOGRAPHICLIB_GEODESIC_ORDER 6
#define nC1   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC1p  GEOGRAPHICLIB_GEODESIC_ORDER
#define nC2   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3   GEOGRAPHICLIB_GEODESIC_ORDER
#define nA3x  nA3
#define nC3   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC3x  ((nC3 * (nC3 - 1)) / 2)
#define nC4   GEOGRAPHICLIB_GEODESIC_ORDER
#define nC4x  ((nC4 * (nC4 + 1)) / 2)

typedef double real;
typedef int boolx;

static unsigned init = 0;
static const int FALSE = 0;
static const int TRUE = 1;
static unsigned digits, maxit1, maxit2;
static real epsilon, realmin, pi, degree, NaN,
  tiny, tol0, tol1, tol2, tolb, xthresh;

static void Init() {
  if (!init) {
#if defined(__DBL_MANT_DIG__)
    digits = __DBL_MANT_DIG__;
#else
    digits = 53;
#endif
#if defined(__DBL_EPSILON__)
    epsilon = __DBL_EPSILON__;
#else
    epsilon = pow(0.5, digits - 1);
#endif
#if defined(__DBL_MIN__)
    realmin = __DBL_MIN__;
#else
    realmin = pow(0.5, 1022);
#endif
#if defined(M_PI)
    pi = M_PI;
#else
    pi = atan2(0.0, -1.0);
#endif
    maxit1 = 20;
    maxit2 = maxit1 + digits + 10;
    tiny = sqrt(realmin);
    tol0 = epsilon;
    /* Increase multiplier in defn of tol1 from 100 to 200 to fix inverse case
     * 52.784459512564 0 -52.784459512563990912 179.634407464943777557
     * which otherwise failed for Visual Studio 10 (Release and Debug) */
    tol1 = 200 * tol0;
    tol2 = sqrt(tol0);
    /* Check on bisection interval */
    tolb = tol0 * tol2;
    xthresh = 1000 * tol2;
    degree = pi/180;
    NaN = sqrt(-1.0);
    init = 1;
  }
}

enum captype {
  CAP_NONE = 0U,
  CAP_C1   = 1U<<0,
  CAP_C1p  = 1U<<1,
  CAP_C2   = 1U<<2,
  CAP_C3   = 1U<<3,
  CAP_C4   = 1U<<4,
  CAP_ALL  = 0x1FU,
  OUT_ALL  = 0x7F80U
};

static real sq(real x) { return x * x; }
static real log1px(real x) {
  volatile real
    y = 1 + x,
    z = y - 1;
  /* Here's the explanation for this magic: y = 1 + z, exactly, and z
   * approx x, thus log(y)/z (which is nearly constant near z = 0) returns
   * a good approximation to the true log(1 + x)/x.  The multiplication x *
   * (log(y)/z) introduces little additional error. */
  return z == 0 ? x : x * log(y) / z;
}

static real atanhx(real x) {
  real y = fabs(x);             /* Enforce odd parity */
  y = log1px(2 * y/(1 - y))/2;
  return x < 0 ? -y : y;
}

static real hypotx(real x, real y)
{ return sqrt(x * x + y * y); }

static real cbrtx(real x) {
  real y = pow(fabs(x), 1/(real)(3)); /* Return the real cube root */
  return x < 0 ? -y : y;
}

static real minx(real x, real y)
{ return x < y ? x : y; }

static real maxx(real x, real y)
{ return x > y ? x : y; }

static void swapx(real* x, real* y)
{ real t = *x; *x = *y; *y = t; }

static void SinCosNorm(real* sinx, real* cosx) {
  real r = hypotx(*sinx, *cosx);
  *sinx /= r;
  *cosx /= r;
}

static real AngNormalize(real x)
{ return x >= 180 ? x - 360 : (x < -180 ? x + 360 : x); }
static real AngNormalize2(real x)
{ return AngNormalize(fmod(x, (real)(360))); }

static real AngRound(real x) {
  const real z = (real)(0.0625); /* 1/16 */
  volatile real y = fabs(x);
  /* The compiler mustn't "simplify" z - (z - y) to y */
  y = y < z ? z - (z - y) : y;
  return x < 0 ? -y : y;
}

static void A3coeff(struct Geodesic* g);
static void C3coeff(struct Geodesic* g);
static void C4coeff(struct Geodesic* g);
static real SinCosSeries(boolx sinp,
                         real sinx, real cosx,
                         const real c[], int n);
static void Lengths(const struct Geodesic* g,
                    real eps, real sig12,
                    real ssig1, real csig1, real dn1,
                    real ssig2, real csig2, real dn2,
                    real cbet1, real cbet2,
                    real* ps12b, real* pm12b, real* pm0,
                    boolx scalep, real* pM12, real* pM21,
                    /* Scratch areas of the right size */
                    real C1a[], real C2a[]);
static real Astroid(real x, real y);
static real InverseStart(const struct Geodesic* g,
                         real sbet1, real cbet1, real dn1,
                         real sbet2, real cbet2, real dn2,
                         real lam12,
                         real* psalp1, real* pcalp1,
                         /* Only updated if return val >= 0 */
                         real* psalp2, real* pcalp2,
                         /* Scratch areas of the right size */
                         real C1a[], real C2a[]);
static real Lambda12(const struct Geodesic* g,
                     real sbet1, real cbet1, real dn1,
                     real sbet2, real cbet2, real dn2,
                     real salp1, real calp1,
                     real* psalp2, real* pcalp2,
                     real* psig12,
                     real* pssig1, real* pcsig1,
                     real* pssig2, real* pcsig2,
                     real* peps, real* pdomg12,
                     boolx diffp, real* pdlam12,
                     /* Scratch areas of the right size */
                     real C1a[], real C2a[], real C3a[]);
static real A3f(const struct Geodesic* g, real eps);
static void C3f(const struct Geodesic* g, real eps, real c[]);
static void C4f(const struct Geodesic* g, real eps, real c[]);
static real A1m1f(real eps);
static void C1f(real eps, real c[]);
static void C1pf(real eps, real c[]);
static real A2m1f(real eps);
static void C2f(real eps, real c[]);

void GeodesicInit(struct Geodesic* g, real a, real f) {
  if (!init) Init();
  g->a = a;
  g->f = f <= 1 ? f : 1/f;
  g->f1 = 1 - g->f;
  g->e2 = g->f * (2 - g->f);
  g->ep2 = g->e2 / sq(g->f1);   /* e2 / (1 - e2) */
  g->n = g->f / ( 2 - g->f);
  g->b = g->a * g->f1;
  g->c2 = (sq(g->a) + sq(g->b) *
           (g->e2 == 0 ? 1 :
            (g->e2 > 0 ? atanhx(sqrt(g->e2)) : atan(sqrt(-g->e2))) /
            sqrt(fabs(g->e2))))/2; /* authalic radius squared */
  /* The sig12 threshold for "really short" */
  g->etol2 = 0.01 * tol2 / maxx((real)(0.1), sqrt(fabs(g->e2)));
  A3coeff(g);
  C3coeff(g);
  C4coeff(g);
}

void GeodesicLineInit(struct GeodesicLine* l,
                      const struct Geodesic* g,
                      real lat1, real lon1, real azi1, unsigned caps) {
  real alp1, cbet1, sbet1, phi, eps;
  l->a = g->a;
  l->f = g->f;
  l->b = g->b;
  l->c2 = g->c2;
  l->f1 = g->f1;
  /* If caps is 0 assume the standard direct calculation */
  l->caps = (caps ? caps : DISTANCE_IN | LONGITUDE) |
    LATITUDE | AZIMUTH; /* Always allow latitude and azimuth */

  azi1 = AngNormalize(azi1);
  lon1 = AngNormalize(lon1);
  if (lat1 == 90) {
    lon1 += lon1 < 0 ? 180 : -180;
    lon1 = AngNormalize(lon1 - azi1);
    azi1 = -180;
  } else if (lat1 == -90) {
    lon1 = AngNormalize(lon1 + azi1);
    azi1 = 0;
  } else {
    /* Guard against underflow in salp0 */
    azi1 = AngRound(azi1);
  }

  l->lat1 = lat1;
  l->lon1 = lon1;
  l->azi1 = azi1;
  /* alp1 is in [0, pi] */
  alp1 = azi1 * degree;
  /* Enforce sin(pi) == 0 and cos(pi/2) == 0.  Better to face the ensuing
   * problems directly than to skirt them. */
  l->salp1 =     azi1  == -180 ? 0 : sin(alp1);
  l->calp1 = fabs(azi1) ==   90 ? 0 : cos(alp1);
  phi = lat1 * degree;
  /* Ensure cbet1 = +epsilon at poles */
  sbet1 = l->f1 * sin(phi);
  cbet1 = fabs(lat1) == 90 ? tiny : cos(phi);
  SinCosNorm(&sbet1, &cbet1);
  l->dn1 = sqrt(1 + g->ep2 * sq(sbet1));

  /* Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0), */
  l->salp0 = l->salp1 * cbet1; /* alp0 in [0, pi/2 - |bet1|] */
  /* Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
   * is slightly better (consider the case salp1 = 0). */
  l->calp0 = hypotx(l->calp1, l->salp1 * sbet1);
  /* Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
   * sig = 0 is nearest northward crossing of equator.
   * With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
   * With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
   * With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
   * Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
   * With alp0 in (0, pi/2], quadrants for sig and omg coincide.
   * No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
   * With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi. */
  l->ssig1 = sbet1; l->somg1 = l->salp0 * sbet1;
  l->csig1 = l->comg1 = sbet1 != 0 || l->calp1 != 0 ? cbet1 * l->calp1 : 1;
  SinCosNorm(&l->ssig1, &l->csig1); /* sig1 in (-pi, pi] */
  /* SinCosNorm(somg1, comg1); -- don't need to normalize! */

  l->k2 = sq(l->calp0) * g->ep2;
  eps = l->k2 / (2 * (1 + sqrt(1 + l->k2)) + l->k2);

  if (l->caps & CAP_C1) {
    real s, c;
    l->A1m1 = A1m1f(eps);
    C1f(eps, l->C1a);
    l->B11 = SinCosSeries(TRUE, l->ssig1, l->csig1, l->C1a, nC1);
    s = sin(l->B11); c = cos(l->B11);
    /* tau1 = sig1 + B11 */
    l->stau1 = l->ssig1 * c + l->csig1 * s;
    l->ctau1 = l->csig1 * c - l->ssig1 * s;
    /* Not necessary because C1pa reverts C1a
     *    B11 = -SinCosSeries(TRUE, stau1, ctau1, C1pa, nC1p); */
  }

  if (l->caps & CAP_C1p)
    C1pf(eps, l->C1pa);

  if (l->caps & CAP_C2) {
    l->A2m1 = A2m1f(eps);
    C2f(eps, l->C2a);
    l->B21 = SinCosSeries(TRUE, l->ssig1, l->csig1, l->C2a, nC2);
  }

  if (l->caps & CAP_C3) {
    C3f(g, eps, l->C3a);
    l->A3c = -l->f * l->salp0 * A3f(g, eps);
    l->B31 = SinCosSeries(TRUE, l->ssig1, l->csig1, l->C3a, nC3-1);
  }

  if (l->caps & CAP_C4) {
    C4f(g, eps, l->C4a);
    /* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0) */
    l->A4 = sq(l->a) * l->calp0 * l->salp0 * g->e2;
    l->B41 = SinCosSeries(FALSE, l->ssig1, l->csig1, l->C4a, nC4);
  }
}

real GenPosition(const struct GeodesicLine* l,
                 boolx arcmode, real s12_a12,
                 real* plat2, real* plon2, real* pazi2,
                 real* ps12, real* pm12,
                 real* pM12, real* pM21,
                 real* pS12) {
  real lat2 = 0, lon2 = 0, azi2 = 0, s12 = 0,
    m12 = 0, M12 = 0, M21 = 0, S12 = 0;
  /* Avoid warning about uninitialized B12. */
  real sig12, ssig12, csig12, B12 = 0, AB1 = 0;
  real omg12, lam12, lon12;
  real ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2, dn2;
  unsigned outmask =
    (plat2 ? LATITUDE : 0U) |
    (plon2 ? LONGITUDE : 0U) |
    (pazi2 ? AZIMUTH : 0U) |
    (ps12 ? DISTANCE : 0U) |
    (pm12 ? REDUCEDLENGTH : 0U) |
    (pM12 || pM21 ? GEODESICSCALE : 0U) |
    (pS12 ? AREA : 0U);

  outmask &= l->caps & OUT_ALL;
  if (!( TRUE /*Init()*/ && (arcmode || (l->caps & DISTANCE_IN & OUT_ALL)) ))
    /* Uninitialized or impossible distance calculation requested */
    return NaN;

  if (arcmode) {
    real s12a;
    /* Interpret s12_a12 as spherical arc length */
    sig12 = s12_a12 * degree;
    s12a = fabs(s12_a12);
    s12a -= 180 * floor(s12a / 180);
    ssig12 = s12a ==  0 ? 0 : sin(sig12);
    csig12 = s12a == 90 ? 0 : cos(sig12);
  } else {
    /* Interpret s12_a12 as distance */
    real
      tau12 = s12_a12 / (l->b * (1 + l->A1m1)),
      s = sin(tau12),
      c = cos(tau12);
    /* tau2 = tau1 + tau12 */
    B12 = - SinCosSeries(TRUE,
                         l->stau1 * c + l->ctau1 * s,
                         l->ctau1 * c - l->stau1 * s,
                         l->C1pa, nC1p);
    sig12 = tau12 - (B12 - l->B11);
    ssig12 = sin(sig12); csig12 = cos(sig12);
    if (fabs(l->f) > 0.01) {
      /* Reverted distance series is inaccurate for |f| > 1/100, so correct
       * sig12 with 1 Newton iteration.  The following table shows the
       * approximate maximum error for a = WGS_a() and various f relative to
       * GeodesicExact.
       *     erri = the error in the inverse solution (nm)
       *     errd = the error in the direct solution (series only) (nm)
       *     errda = the error in the direct solution (series + 1 Newton) (nm)
       *
       *       f     erri  errd errda
       *     -1/5    12e6 1.2e9  69e6
       *     -1/10  123e3  12e6 765e3
       *     -1/20   1110 108e3  7155
       *     -1/50  18.63 200.9 27.12
       *     -1/100 18.63 23.78 23.37
       *     -1/150 18.63 21.05 20.26
       *      1/150 22.35 24.73 25.83
       *      1/100 22.35 25.03 25.31
       *      1/50  29.80 231.9 30.44
       *      1/20   5376 146e3  10e3
       *      1/10  829e3  22e6 1.5e6
       *      1/5   157e6 3.8e9 280e6 */
      real
        ssig2 = l->ssig1 * csig12 + l->csig1 * ssig12,
        csig2 = l->csig1 * csig12 - l->ssig1 * ssig12,
        serr;
      B12 = SinCosSeries(TRUE, ssig2, csig2, l->C1a, nC1);
      serr = (1 + l->A1m1) * (sig12 + (B12 - l->B11)) - s12_a12 / l->b;
      sig12 = sig12 - serr / sqrt(1 + l->k2 * sq(ssig2));
      ssig12 = sin(sig12); csig12 = cos(sig12);
      /* Update B12 below */
    }
  }

  /* sig2 = sig1 + sig12 */
  ssig2 = l->ssig1 * csig12 + l->csig1 * ssig12;
  csig2 = l->csig1 * csig12 - l->ssig1 * ssig12;
  dn2 = sqrt(1 + l->k2 * sq(ssig2));
  if (outmask & (DISTANCE | REDUCEDLENGTH | GEODESICSCALE)) {
    if (arcmode || fabs(l->f) > 0.01)
      B12 = SinCosSeries(TRUE, ssig2, csig2, l->C1a, nC1);
    AB1 = (1 + l->A1m1) * (B12 - l->B11);
  }
  /* sin(bet2) = cos(alp0) * sin(sig2) */
  sbet2 = l->calp0 * ssig2;
  /* Alt: cbet2 = hypot(csig2, salp0 * ssig2); */
  cbet2 = hypotx(l->salp0, l->calp0 * csig2);
  if (cbet2 == 0)
    /* I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case */
    cbet2 = csig2 = tiny;
  /* tan(omg2) = sin(alp0) * tan(sig2) */
  somg2 = l->salp0 * ssig2; comg2 = csig2;  /* No need to normalize */
  /* tan(alp0) = cos(sig2)*tan(alp2) */
  salp2 = l->salp0; calp2 = l->calp0 * csig2; /* No need to normalize */
  /* omg12 = omg2 - omg1 */
  omg12 = atan2(somg2 * l->comg1 - comg2 * l->somg1,
                comg2 * l->comg1 + somg2 * l->somg1);

  if (outmask & DISTANCE)
    s12 = arcmode ? l->b * ((1 + l->A1m1) * sig12 + AB1) : s12_a12;

  if (outmask & LONGITUDE) {
    lam12 = omg12 + l->A3c *
      ( sig12 + (SinCosSeries(TRUE, ssig2, csig2, l->C3a, nC3-1)
                 - l->B31));
    lon12 = lam12 / degree;
    /* Use AngNormalize2 because longitude might have wrapped multiple
     * times. */
    lon12 = AngNormalize2(lon12);
    lon2 = AngNormalize(l->lon1 + lon12);
  }

  if (outmask & LATITUDE)
    lat2 = atan2(sbet2, l->f1 * cbet2) / degree;

  if (outmask & AZIMUTH)
    /* minus signs give range [-180, 180). 0- converts -0 to +0. */
    azi2 = 0 - atan2(-salp2, calp2) / degree;

  if (outmask & (REDUCEDLENGTH | GEODESICSCALE)) {
    real
      B22 = SinCosSeries(TRUE, ssig2, csig2, l->C2a, nC2),
      AB2 = (1 + l->A2m1) * (B22 - l->B21),
      J12 = (l->A1m1 - l->A2m1) * sig12 + (AB1 - AB2);
    if (outmask & REDUCEDLENGTH)
      /* Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
       * accurate cancellation in the case of coincident points. */
      m12 = l->b * ((dn2 * (l->csig1 * ssig2) - l->dn1 * (l->ssig1 * csig2))
                    - l->csig1 * csig2 * J12);
    if (outmask & GEODESICSCALE) {
      real t = l->k2 * (ssig2 - l->ssig1) * (ssig2 + l->ssig1) / (l->dn1 + dn2);
      M12 = csig12 + (t *  ssig2 -  csig2 * J12) * l->ssig1 / l->dn1;
      M21 = csig12 - (t * l->ssig1 - l->csig1 * J12) *  ssig2 /  dn2;
    }
  }

  if (outmask & AREA) {
    real
      B42 = SinCosSeries(FALSE, ssig2, csig2, l->C4a, nC4);
    real salp12, calp12;
    if (l->calp0 == 0 || l->salp0 == 0) {
      /* alp12 = alp2 - alp1, used in atan2 so no need to normalized */
      salp12 = salp2 * l->calp1 - calp2 * l->salp1;
      calp12 = calp2 * l->calp1 + salp2 * l->salp1;
      /* The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
       * salp12 = -0 and alp12 = -180.  However this depends on the sign being
       * attached to 0 correctly.  The following ensures the correct
       * behavior. */
      if (salp12 == 0 && calp12 < 0) {
        salp12 = tiny * l->calp1;
        calp12 = -1;
      }
    } else {
      /* tan(alp) = tan(alp0) * sec(sig)
       * tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
       * = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
       * If csig12 > 0, write
       *   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
       * else
       *   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
       * No need to normalize */
      salp12 = l->calp0 * l->salp0 *
        (csig12 <= 0 ? l->csig1 * (1 - csig12) + ssig12 * l->ssig1 :
         ssig12 * (l->csig1 * ssig12 / (1 + csig12) + l->ssig1));
      calp12 = sq(l->salp0) + sq(l->calp0) * l->csig1 * csig2;
    }
    S12 = l->c2 * atan2(salp12, calp12) + l->A4 * (B42 - l->B41);
  }

  if (outmask & LATITUDE)
    *plat2 = lat2;
  if (outmask & LONGITUDE)
    *plon2 = lon2;
  if (outmask & AZIMUTH)
    *pazi2 = azi2;
  if (outmask & DISTANCE)
    *ps12 = s12;
  if (outmask & REDUCEDLENGTH)
    *pm12 = m12;
  if (outmask & GEODESICSCALE) {
    if (pM12) *pM12 = M12;
    if (pM21) *pM21 = M21;
  }
  if (outmask & AREA)
    *pS12 = S12;

  return arcmode ? s12_a12 : sig12 / degree;
}

void Position(const struct GeodesicLine* l, real s12,
              real* plat2, real* plon2, real* pazi2) {
  GenPosition(l, FALSE, s12, plat2, plon2, pazi2, 0, 0, 0, 0, 0);
}

real GenDirect(const struct Geodesic* g,
               real lat1, real lon1, real azi1,
               boolx arcmode, real s12_a12,
               real* plat2, real* plon2, real* pazi2,
               real* ps12, real* pm12, real* pM12, real* pM21,
               real* pS12) {
  struct GeodesicLine l;
  unsigned outmask =
    (plat2 ? LATITUDE : 0U) |
    (plon2 ? LONGITUDE : 0U) |
    (pazi2 ? AZIMUTH : 0U) |
    (ps12 ? DISTANCE : 0U) |
    (pm12 ? REDUCEDLENGTH : 0U) |
    (pM12 || pM21 ? GEODESICSCALE : 0U) |
    (pS12 ? AREA : 0U);

  GeodesicLineInit(&l, g, lat1, lon1, azi1,
                   /* Automatically supply DISTANCE_IN if necessary */
                   outmask | (arcmode ? NONE : DISTANCE_IN));
  return GenPosition(&l, arcmode, s12_a12,
                     plat2, plon2, pazi2, ps12, pm12, pM12, pM21, pS12);
}

void Direct(const struct Geodesic* g,
            real lat1, real lon1, real azi1,
            real s12,
            real* plat2, real* plon2, real* pazi2) {
  GenDirect(g, lat1, lon1, azi1, FALSE, s12, plat2, plon2, pazi2,
            0, 0, 0, 0, 0);
}

real GenInverse(const struct Geodesic* g,
                real lat1, real lon1, real lat2, real lon2,
                real* ps12, real* pazi1, real* pazi2,
                real* pm12, real* pM12, real* pM21, real* pS12) {
  real s12 = 0, azi1 = 0, azi2 = 0, m12 = 0, M12 = 0, M21 = 0, S12 = 0;
  real lon12;
  int latsign, lonsign, swapp;
  real phi, sbet1, cbet1, sbet2, cbet2, s12x = 0, m12x = 0;
  real dn1, dn2, lam12, slam12, clam12;
  real a12 = 0, sig12, calp1 = 0, salp1 = 0, calp2 = 0, salp2 = 0;
  /* index zero elements of these arrays are unused */
  real C1a[nC1 + 1], C2a[nC2 + 1], C3a[nC3];
  boolx meridian;
  real omg12 = 0;

  unsigned outmask =
    (ps12 ? DISTANCE : 0U) |
    (pazi1 || pazi2 ? AZIMUTH : 0U) |
    (pm12 ? REDUCEDLENGTH : 0U) |
    (pM12 || pM21 ? GEODESICSCALE : 0U) |
    (pS12 ? AREA : 0U);

  outmask &= OUT_ALL;
  lon12 = AngNormalize(AngNormalize(lon2) -
                            AngNormalize(lon1));
  /* If very close to being on the same meridian, then make it so.
   * Not sure this is necessary... */
  lon12 = AngRound(lon12);
  /* Make longitude difference positive. */
  lonsign = lon12 >= 0 ? 1 : -1;
  lon12 *= lonsign;
  if (lon12 == 180)
    lonsign = 1;
  /* If really close to the equator, treat as on equator. */
  lat1 = AngRound(lat1);
  lat2 = AngRound(lat2);
  /* Swap points so that point with higher (abs) latitude is point 1 */
  swapp = fabs(lat1) >= fabs(lat2) ? 1 : -1;
  if (swapp < 0) {
    lonsign *= -1;
    swapx(&lat1, &lat2);
  }
  /* Make lat1 <= 0 */
  latsign = lat1 < 0 ? 1 : -1;
  lat1 *= latsign;
  lat2 *= latsign;
  /* Now we have
   *
   *     0 <= lon12 <= 180
   *     -90 <= lat1 <= 0
   *     lat1 <= lat2 <= -lat1
   *
   * longsign, swapp, latsign register the transformation to bring the
   * coordinates to this canonical form.  In all cases, 1 means no change was
   * made.  We make these transformations so that there are few cases to
   * check, e.g., on verifying quadrants in atan2.  In addition, this
   * enforces some symmetries in the results returned. */

  phi = lat1 * degree;
  /* Ensure cbet1 = +epsilon at poles */
  sbet1 = g->f1 * sin(phi);
  cbet1 = lat1 == -90 ? tiny : cos(phi);
  SinCosNorm(&sbet1, &cbet1);

  phi = lat2 * degree;
  /* Ensure cbet2 = +epsilon at poles */
  sbet2 = g->f1 * sin(phi);
  cbet2 = fabs(lat2) == 90 ? tiny : cos(phi);
  SinCosNorm(&sbet2, &cbet2);

  /* If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
   * |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
   * a better measure.  This logic is used in assigning calp2 in Lambda12.
   * Sometimes these quantities vanish and in that case we force bet2 = +/-
   * bet1 exactly.  An example where is is necessary is the inverse problem
   * 48.522876735459 0 -48.52287673545898293 179.599720456223079643
   * which failed with Visual Studio 10 (Release and Debug) */

  if (cbet1 < -sbet1) {
    if (cbet2 == cbet1)
      sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
  } else {
    if (fabs(sbet2) == -sbet1)
      cbet2 = cbet1;
  }

  dn1 = sqrt(1 + g->ep2 * sq(sbet1));
  dn2 = sqrt(1 + g->ep2 * sq(sbet2));

  lam12 = lon12 * degree;
  slam12 = lon12 == 180 ? 0 : sin(lam12);
  clam12 = cos(lam12);      /* lon12 == 90 isn't interesting */


  meridian = lat1 == -90 || slam12 == 0;

  if (meridian) {

    /* Endpoints are on a single full meridian, so the geodesic might lie on
     * a meridian. */

    real ssig1, csig1, ssig2, csig2;
    calp1 = clam12; salp1 = slam12; /* Head to the target longitude */
    calp2 = 1; salp2 = 0;           /* At the target we're heading north */

    /* tan(bet) = tan(sig) * cos(alp) */
    ssig1 = sbet1; csig1 = calp1 * cbet1;
    ssig2 = sbet2; csig2 = calp2 * cbet2;

    /* sig12 = sig2 - sig1 */
    sig12 = atan2(maxx(csig1 * ssig2 - ssig1 * csig2, (real)(0)),
                  csig1 * csig2 + ssig1 * ssig2);
    {
      real dummy;
      Lengths(g, g->n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
              cbet1, cbet2, &s12x, &m12x, &dummy,
              (outmask & GEODESICSCALE) != 0U, &M12, &M21, C1a, C2a);
    }
    /* Add the check for sig12 since zero length geodesics might yield m12 <
     * 0.  Test case was
     *
     *    echo 20.001 0 20.001 0 | Geod -i
     *
     * In fact, we will have sig12 > pi/2 for meridional geodesic which is
     * not a shortest path. */
    if (sig12 < 1 || m12x >= 0) {
      m12x *= g->b;
      s12x *= g->b;
      a12 = sig12 / degree;
    } else
      /* m12 < 0, i.e., prolate and too close to anti-podal */
      meridian = FALSE;
  }

  if (!meridian &&
      sbet1 == 0 &&           /* and sbet2 == 0 */
      /* Mimic the way Lambda12 works with calp1 = 0 */
      (g->f <= 0 || lam12 <= pi - g->f * pi)) {

    /* Geodesic runs along equator */
    calp1 = calp2 = 0; salp1 = salp2 = 1;
    s12x = g->a * lam12;
    sig12 = omg12 = lam12 / g->f1;
    m12x = g->b * sin(sig12);
    if (outmask & GEODESICSCALE)
      M12 = M21 = cos(sig12);
    a12 = lon12 / g->f1;

  } else if (!meridian) {

    /* Now point1 and point2 belong within a hemisphere bounded by a
     * meridian and geodesic is neither meridional or equatorial. */

    /* Figure a starting point for Newton's method */
    sig12 = InverseStart(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
                         lam12,
                         &salp1, &calp1, &salp2, &calp2,
                         C1a, C2a);

    if (sig12 >= 0) {
      /* Short lines (InverseStart sets salp2, calp2) */
      real dnm = (dn1 + dn2) / 2;
      s12x = sig12 * g->b * dnm;
      m12x = sq(dnm) * g->b * sin(sig12 / dnm);
      if (outmask & GEODESICSCALE)
        M12 = M21 = cos(sig12 / dnm);
      a12 = sig12 / degree;
      omg12 = lam12 / (g->f1 * dnm);
    } else {

      /* Newton's method.  This is a straightforward solution of f(alp1) =
       * lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
       * root in the interval (0, pi) and its derivative is positive at the
       * root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
       * alp1.  During the course of the iteration, a range (alp1a, alp1b) is
       * maintained which brackets the root and with each evaluation of
       * f(alp) the range is shrunk, if possible.  Newton's method is
       * restarted whenever the derivative of f is negative (because the new
       * value of alp1 is then further from the solution) or if the new
       * estimate of alp1 lies outside (0,pi); in this case, the new starting
       * guess is taken to be (alp1a + alp1b) / 2. */
      real ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0;
      unsigned numit = 0;
      /* Bracketing range */
      real salp1a = tiny, calp1a = 1, salp1b = tiny, calp1b = -1;
      boolx tripn, tripb;
      for (tripn = FALSE, tripb = FALSE; numit < maxit2; ++numit) {
        /* the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
         * WGS84 and random input: mean = 2.85, sd = 0.60 */
        real dv,
          v = (Lambda12(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
                        &salp2, &calp2, &sig12, &ssig1, &csig1, &ssig2, &csig2,
                        &eps, &omg12, numit < maxit1, &dv, C1a, C2a, C3a)
               - lam12);
        /* 2 * tol0 is approximately 1 ulp for a number in [0, pi]. */
        if (tripb || fabs(v) < (tripn ? 8 : 2) * tol0) break;
        /* Update bracketing values */
        if (v > 0 && (numit > maxit1 || calp1/salp1 > calp1b/salp1b)) {
          salp1b = salp1; calp1b = calp1;
        } else if (numit > maxit1 || calp1/salp1 < calp1a/salp1a) {
          salp1a = salp1; calp1a = calp1;
        }
        if (numit < maxit1 && dv > 0) {
          real
            dalp1 = -v/dv;
          real
            sdalp1 = sin(dalp1), cdalp1 = cos(dalp1),
            nsalp1 = salp1 * cdalp1 + calp1 * sdalp1;
          if (nsalp1 > 0 && fabs(dalp1) < pi) {
            calp1 = calp1 * cdalp1 - salp1 * sdalp1;
            salp1 = nsalp1;
            SinCosNorm(&salp1, &calp1);
            /* In some regimes we don't get quadratic convergence because
             * slope -> 0.  So use convergence conditions based on epsilon
             * instead of sqrt(epsilon). */
            tripn = fabs(v) <= 16 * tol0;
            continue;
          }
        }
        /* Either dv was not postive or updated value was outside legal
         * range.  Use the midpoint of the bracket as the next estimate.
         * This mechanism is not needed for the WGS84 ellipsoid, but it does
         * catch problems with more eccentric ellipsoids.  Its efficacy is
         * such for the WGS84 test set with the starting guess set to alp1 =
         * 90deg:
         * the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
         * WGS84 and random input: mean = 4.74, sd = 0.99 */
        salp1 = (salp1a + salp1b)/2;
        calp1 = (calp1a + calp1b)/2;
        SinCosNorm(&salp1, &calp1);
        tripn = FALSE;
        tripb = (fabs(salp1a - salp1) + (calp1a - calp1) < tolb ||
                 fabs(salp1 - salp1b) + (calp1 - calp1b) < tolb);
      }
      {
        real dummy;
        Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
                cbet1, cbet2, &s12x, &m12x, &dummy,
                (outmask & GEODESICSCALE) != 0U, &M12, &M21, C1a, C2a);
      }
      m12x *= g->b;
      s12x *= g->b;
      a12 = sig12 / degree;
      omg12 = lam12 - omg12;
    }
  }

  if (outmask & DISTANCE)
    s12 = 0 + s12x;             /* Convert -0 to 0 */

  if (outmask & REDUCEDLENGTH)
    m12 = 0 + m12x;             /* Convert -0 to 0 */

  if (outmask & AREA) {
    real
      /* From Lambda12: sin(alp1) * cos(bet1) = sin(alp0) */
      salp0 = salp1 * cbet1,
      calp0 = hypotx(calp1, salp1 * sbet1); /* calp0 > 0 */
    real alp12;
    if (calp0 != 0 && salp0 != 0) {
      real
        /* From Lambda12: tan(bet) = tan(sig) * cos(alp) */
        ssig1 = sbet1, csig1 = calp1 * cbet1,
        ssig2 = sbet2, csig2 = calp2 * cbet2,
        k2 = sq(calp0) * g->ep2,
        eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2),
        /* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0). */
        A4 = sq(g->a) * calp0 * salp0 * g->e2;
      real C4a[nC4];
      real B41, B42;
      SinCosNorm(&ssig1, &csig1);
      SinCosNorm(&ssig2, &csig2);
      C4f(g, eps, C4a);
      B41 = SinCosSeries(FALSE, ssig1, csig1, C4a, nC4);
      B42 = SinCosSeries(FALSE, ssig2, csig2, C4a, nC4);
      S12 = A4 * (B42 - B41);
    } else
      /* Avoid problems with indeterminate sig1, sig2 on equator */
      S12 = 0;

    if (!meridian &&
        omg12 < (real)(0.75) * pi &&   /* Long difference too big */
        sbet2 - sbet1 < (real)(1.75)) { /* Lat difference too big */
      /* Use tan(Gamma/2) = tan(omg12/2)
       * * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
       * with tan(x/2) = sin(x)/(1+cos(x)) */
      real
        somg12 = sin(omg12), domg12 = 1 + cos(omg12),
        dbet1 = 1 + cbet1, dbet2 = 1 + cbet2;
      alp12 = 2 * atan2( somg12 * ( sbet1 * dbet2 + sbet2 * dbet1 ),
                         domg12 * ( sbet1 * sbet2 + dbet1 * dbet2 ) );
    } else {
      /* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
      real
        salp12 = salp2 * calp1 - calp2 * salp1,
        calp12 = calp2 * calp1 + salp2 * salp1;
      /* The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
       * salp12 = -0 and alp12 = -180.  However this depends on the sign
       * being attached to 0 correctly.  The following ensures the correct
       * behavior. */
      if (salp12 == 0 && calp12 < 0) {
        salp12 = tiny * calp1;
        calp12 = -1;
      }
      alp12 = atan2(salp12, calp12);
    }
    S12 += g->c2 * alp12;
    S12 *= swapp * lonsign * latsign;
    /* Convert -0 to 0 */
    S12 += 0;
  }

  /* Convert calp, salp to azimuth accounting for lonsign, swapp, latsign. */
  if (swapp < 0) {
    swapx(&salp1, &salp2);
    swapx(&calp1, &calp2);
    if (outmask & GEODESICSCALE)
      swapx(&M12, &M21);
  }

  salp1 *= swapp * lonsign; calp1 *= swapp * latsign;
  salp2 *= swapp * lonsign; calp2 *= swapp * latsign;

  if (outmask & AZIMUTH) {
    /* minus signs give range [-180, 180). 0- converts -0 to +0. */
    azi1 = 0 - atan2(-salp1, calp1) / degree;
    azi2 = 0 - atan2(-salp2, calp2) / degree;
  }

  if (outmask & DISTANCE)
    *ps12 = s12;
  if (outmask & AZIMUTH) {
    if (pazi1) *pazi1 = azi1;
    if (pazi2) *pazi2 = azi2;
  }
  if (outmask & REDUCEDLENGTH)
    *pm12 = m12;
  if (outmask & GEODESICSCALE) {
    if (pM12) *pM12 = M12;
    if (pM21) *pM21 = M21;
  }
  if (outmask & AREA)
    *pS12 = S12;

  /* Returned value in [0, 180] */
  return a12;
}

void Inverse(const struct Geodesic* g,
             double lat1, double lon1, double lat2, double lon2,
             double* ps12, double* pazi1, double* pazi2) {
  GenInverse(g, lat1, lon1, lat2, lon2, ps12, pazi1, pazi2, 0, 0, 0, 0);
}

real SinCosSeries(boolx sinp,
                  real sinx, real cosx,
                  const real c[], int n) {
  /* Evaluate
   * y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
   *            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
   * using Clenshaw summation.  N.B. c[0] is unused for sin series
   * Approx operation count = (n + 5) mult and (2 * n + 2) add */
  real ar, y0, y1;
  c += (n + sinp);              /* Point to one beyond last element */
  ar = 2 * (cosx - sinx) * (cosx + sinx); /* 2 * cos(2 * x) */
  y0 = n & 1 ? *--c : 0; y1 = 0;          /* accumulators for sum */
  /* Now n is even */
  n /= 2;
  while (n--) {
    /* Unroll loop x 2, so accumulators return to their original role */
    y1 = ar * y0 - y1 + *--c;
    y0 = ar * y1 - y0 + *--c;
  }
  return sinp
    ? 2 * sinx * cosx * y0      /* sin(2 * x) * y0 */
    : cosx * (y0 - y1);         /* cos(x) * (y0 - y1) */
}

void Lengths(const struct Geodesic* g,
             real eps, real sig12,
             real ssig1, real csig1, real dn1,
             real ssig2, real csig2, real dn2,
             real cbet1, real cbet2,
             real* ps12b, real* pm12b, real* pm0,
             boolx scalep, real* pM12, real* pM21,
             /* Scratch areas of the right size */
             real C1a[], real C2a[]) {
  real s12b = 0, m12b = 0, m0 = 0, M12 = 0, M21 = 0;
  real A1m1, AB1, A2m1, AB2, J12;

  /* Return m12b = (reduced length)/b; also calculate s12b = distance/b,
   * and m0 = coefficient of secular term in expression for reduced length. */
  C1f(eps, C1a);
  C2f(eps, C2a);
  A1m1 = A1m1f(eps);
  AB1 = (1 + A1m1) * (SinCosSeries(TRUE, ssig2, csig2, C1a, nC1) -
                      SinCosSeries(TRUE, ssig1, csig1, C1a, nC1));
  A2m1 = A2m1f(eps);
  AB2 = (1 + A2m1) * (SinCosSeries(TRUE, ssig2, csig2, C2a, nC2) -
                      SinCosSeries(TRUE, ssig1, csig1, C2a, nC2));
  m0 = A1m1 - A2m1;
  J12 = m0 * sig12 + (AB1 - AB2);
  /* Missing a factor of b.
   * Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure accurate
   * cancellation in the case of coincident points. */
  m12b = dn2 * (csig1 * ssig2) - dn1 * (ssig1 * csig2) - csig1 * csig2 * J12;
  /* Missing a factor of b */
  s12b = (1 + A1m1) * sig12 + AB1;
  if (scalep) {
    real csig12 = csig1 * csig2 + ssig1 * ssig2;
    real t = g->ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2);
    M12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
    M21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
  }
  *ps12b = s12b;
  *pm12b = m12b;
  *pm0 = m0;
  if (scalep) {
    *pM12 = M12;
    *pM21 = M21;
  }
}

real Astroid(real x, real y) {
  /* Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
   * This solution is adapted from Geocentric::Reverse. */
  real k;
  real
    p = sq(x),
    q = sq(y),
    r = (p + q - 1) / 6;
  if ( !(q == 0 && r <= 0) ) {
    real
      /* Avoid possible division by zero when r = 0 by multiplying equations
       * for s and t by r^3 and r, resp. */
      S = p * q / 4,            /* S = r^3 * s */
      r2 = sq(r),
      r3 = r * r2,
      /* The discrimant of the quadratic equation for T3.  This is zero on
       * the evolute curve p^(1/3)+q^(1/3) = 1 */
      disc = S * (S + 2 * r3);
    real u = r;
    real v, uv, w;
    if (disc >= 0) {
      real T3 = S + r3, T;
      /* Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
       * of precision due to cancellation.  The result is unchanged because
       * of the way the T is used in definition of u. */
      T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc); /* T3 = (r * t)^3 */
      /* N.B. cbrtx always returns the real root.  cbrtx(-8) = -2. */
      T = cbrtx(T3);            /* T = r * t */
      /* T can be zero; but then r2 / T -> 0. */
      u += T + (T != 0 ? r2 / T : 0);
    } else {
      /* T is complex, but the way u is defined the result is real. */
      real ang = atan2(sqrt(-disc), -(S + r3));
      /* There are three possible cube roots.  We choose the root which
       * avoids cancellation.  Note that disc < 0 implies that r < 0. */
      u += 2 * r * cos(ang / 3);
    }
    v = sqrt(sq(u) + q);              /* guaranteed positive */
    /* Avoid loss of accuracy when u < 0. */
    uv = u < 0 ? q / (v - u) : u + v; /* u+v, guaranteed positive */
    w = (uv - q) / (2 * v);           /* positive? */
    /* Rearrange expression for k to avoid loss of accuracy due to
     * subtraction.  Division by 0 not possible because uv > 0, w >= 0. */
    k = uv / (sqrt(uv + sq(w)) + w);   /* guaranteed positive */
  } else {               /* q == 0 && r <= 0 */
    /* y = 0 with |x| <= 1.  Handle this case directly.
     * for y small, positive root is k = abs(y)/sqrt(1-x^2) */
    k = 0;
  }
  return k;
}

real InverseStart(const struct Geodesic* g,
                  real sbet1, real cbet1, real dn1,
                  real sbet2, real cbet2, real dn2,
                  real lam12,
                  real* psalp1, real* pcalp1,
                  /* Only updated if return val >= 0 */
                  real* psalp2, real* pcalp2,
                  /* Scratch areas of the right size */
                  real C1a[], real C2a[]) {
  real salp1 = 0, calp1 = 0, salp2 = 0, calp2 = 0;

  /* Return a starting point for Newton's method in salp1 and calp1 (function
   * value is -1).  If Newton's method doesn't need to be used, return also
   * salp2 and calp2 and function value is sig12. */
  real
    sig12 = -1,               /* Return value */
    /* bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0] */
    sbet12 = sbet2 * cbet1 - cbet2 * sbet1,
    cbet12 = cbet2 * cbet1 + sbet2 * sbet1;
#if defined(__GNUC__) && __GNUC__ == 4 &&       \
  (__GNUC_MINOR__ < 6 || defined(__MINGW32__))
  /* Volatile declaration needed to fix inverse cases
   * 88.202499451857 0 -88.202499451857 179.981022032992859592
   * 89.262080389218 0 -89.262080389218 179.992207982775375662
   * 89.333123580033 0 -89.333123580032997687 179.99295812360148422
   * which otherwise fail with g++ 4.4.4 x86 -O3 (Linux)
   * and g++ 4.4.0 (mingw) and g++ 4.6.1 (tdm mingw). */
  real sbet12a;
  {
    volatile real xx1 = sbet2 * cbet1;
    volatile real xx2 = cbet2 * sbet1;
    sbet12a = xx1 + xx2;
  }
#else
  real sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
#endif
  boolx shortline = cbet12 >= 0 && sbet12 < (real)(0.5) &&
    lam12 <= pi / 6;
  real
    omg12 = !shortline ? lam12 : lam12 / (g->f1 * (dn1 + dn2) / 2),
    somg12 = sin(omg12), comg12 = cos(omg12);
  real ssig12, csig12;

  salp1 = cbet2 * somg12;
  calp1 = comg12 >= 0 ?
    sbet12 + cbet2 * sbet1 * sq(somg12) / (1 + comg12) :
    sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);

  ssig12 = hypotx(salp1, calp1);
  csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

  if (shortline && ssig12 < g->etol2) {
    /* really short lines */
    salp2 = cbet1 * somg12;
    calp2 = sbet12 - cbet1 * sbet2 * sq(somg12) / (1 + comg12);
    SinCosNorm(&salp2, &calp2);
    /* Set return value */
    sig12 = atan2(ssig12, csig12);
  } else if (fabs(g->n) > (real)(0.1) || /* No astroid calc if too eccentric */
             csig12 >= 0 ||
             ssig12 >= 6 * fabs(g->n) * pi * sq(cbet1)) {
    /* Nothing to do, zeroth order spherical approximation is OK */
  } else {
    /* Scale lam12 and bet2 to x, y coordinate system where antipodal point
     * is at origin and singular point is at y = 0, x = -1. */
    real y, lamscale, betscale;
    /* Volatile declaration needed to fix inverse case
     * 56.320923501171 0 -56.320923501171 179.664747671772880215
     * which otherwise fails with g++ 4.4.4 x86 -O3 */
    volatile real x;
    if (g->f >= 0) {            /* In fact f == 0 does not get here */
      /* x = dlong, y = dlat */
      {
        real
          k2 = sq(sbet1) * g->ep2,
          eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
        lamscale = g->f * cbet1 * A3f(g, eps) * pi;
      }
      betscale = lamscale * cbet1;

      x = (lam12 - pi) / lamscale;
      y = sbet12a / betscale;
    } else {                    /* f < 0 */
      /* x = dlat, y = dlong */
      real
        cbet12a = cbet2 * cbet1 - sbet2 * sbet1,
        bet12a = atan2(sbet12a, cbet12a);
      real m12b, m0, dummy;
      /* In the case of lon12 = 180, this repeats a calculation made in
       * Inverse. */
      Lengths(g, g->n, pi + bet12a,
              sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
              cbet1, cbet2, &dummy, &m12b, &m0, FALSE,
              &dummy, &dummy, C1a, C2a);
      x = -1 + m12b / (cbet1 * cbet2 * m0 * pi);
      betscale = x < -(real)(0.01) ? sbet12a / x :
        -g->f * sq(cbet1) * pi;
      lamscale = betscale / cbet1;
      y = (lam12 - pi) / lamscale;
    }

    if (y > -tol1 && x > -1 - xthresh) {
      /* strip near cut */
      if (g->f >= 0) {
        salp1 = minx((real)(1), -(real)(x)); calp1 = - sqrt(1 - sq(salp1));
      } else {
        calp1 = maxx((real)(x > -tol1 ? 0 : -1), (real)(x));
        salp1 = sqrt(1 - sq(calp1));
      }
    } else {
      /* Estimate alp1, by solving the astroid problem.
       *
       * Could estimate alpha1 = theta + pi/2, directly, i.e.,
       *   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
       *   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
       *
       * However, it's better to estimate omg12 from astroid and use
       * spherical formula to compute alp1.  This reduces the mean number of
       * Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
       * (min 0 max 5).  The changes in the number of iterations are as
       * follows:
       *
       * change percent
       *    1       5
       *    0      78
       *   -1      16
       *   -2       0.6
       *   -3       0.04
       *   -4       0.002
       *
       * The histogram of iterations is (m = number of iterations estimating
       * alp1 directly, n = number of iterations estimating via omg12, total
       * number of trials = 148605):
       *
       *  iter    m      n
       *    0   148    186
       *    1 13046  13845
       *    2 93315 102225
       *    3 36189  32341
       *    4  5396      7
       *    5   455      1
       *    6    56      0
       *
       * Because omg12 is near pi, estimate work with omg12a = pi - omg12 */
      real k = Astroid(x, y);
      real
        omg12a = lamscale * ( g->f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k );
      somg12 = sin(omg12a); comg12 = -cos(omg12a);
      /* Update spherical estimate of alp1 using omg12 instead of lam12 */
      salp1 = cbet2 * somg12;
      calp1 = sbet12a - cbet2 * sbet1 * sq(somg12) / (1 - comg12);
    }
  }
  if (salp1 > 0)              /* Sanity check on starting guess */
    SinCosNorm(&salp1, &calp1);
  else {
    salp1 = 1; calp1 = 0;
  }

  *psalp1 = salp1;
  *pcalp1 = calp1;
  if (sig12 >= 0) {
    *psalp2 = salp2;
    *pcalp2 = calp2;
  }
  return sig12;
}

real Lambda12(const struct Geodesic* g,
              real sbet1, real cbet1, real dn1,
              real sbet2, real cbet2, real dn2,
              real salp1, real calp1,
              real* psalp2, real* pcalp2,
              real* psig12,
              real* pssig1, real* pcsig1,
              real* pssig2, real* pcsig2,
              real* peps, real* pdomg12,
              boolx diffp, real* pdlam12,
              /* Scratch areas of the right size */
              real C1a[], real C2a[], real C3a[]) {
  real salp2 = 0, calp2 = 0, sig12 = 0,
    ssig1 = 0, csig1 = 0, ssig2 = 0, csig2 = 0, eps = 0, domg12 = 0, dlam12 = 0;
  real salp0, calp0;
  real somg1, comg1, somg2, comg2, omg12, lam12;
  real B312, h0, k2;

  if (sbet1 == 0 && calp1 == 0)
    /* Break degeneracy of equatorial line.  This case has already been
     * handled. */
    calp1 = -tiny;

  /* sin(alp1) * cos(bet1) = sin(alp0) */
  salp0 = salp1 * cbet1;
  calp0 = hypotx(calp1, salp1 * sbet1); /* calp0 > 0 */

  /* tan(bet1) = tan(sig1) * cos(alp1)
   * tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1) */
  ssig1 = sbet1; somg1 = salp0 * sbet1;
  csig1 = comg1 = calp1 * cbet1;
  SinCosNorm(&ssig1, &csig1);
  /* SinCosNorm(&somg1, &comg1); -- don't need to normalize! */

  /* Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
   * about this case, since this can yield singularities in the Newton
   * iteration.
   * sin(alp2) * cos(bet2) = sin(alp0) */
  salp2 = cbet2 != cbet1 ? salp0 / cbet2 : salp1;
  /* calp2 = sqrt(1 - sq(salp2))
   *       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
   * and subst for calp0 and rearrange to give (choose positive sqrt
   * to give alp2 in [0, pi/2]). */
  calp2 = cbet2 != cbet1 || fabs(sbet2) != -sbet1 ?
    sqrt(sq(calp1 * cbet1) +
         (cbet1 < -sbet1 ?
          (cbet2 - cbet1) * (cbet1 + cbet2) :
          (sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
    fabs(calp1);
  /* tan(bet2) = tan(sig2) * cos(alp2)
   * tan(omg2) = sin(alp0) * tan(sig2). */
  ssig2 = sbet2; somg2 = salp0 * sbet2;
  csig2 = comg2 = calp2 * cbet2;
  SinCosNorm(&ssig2, &csig2);
  /* SinCosNorm(&somg2, &comg2); -- don't need to normalize! */

  /* sig12 = sig2 - sig1, limit to [0, pi] */
  sig12 = atan2(maxx(csig1 * ssig2 - ssig1 * csig2, (real)(0)),
                csig1 * csig2 + ssig1 * ssig2);

  /* omg12 = omg2 - omg1, limit to [0, pi] */
  omg12 = atan2(maxx(comg1 * somg2 - somg1 * comg2, (real)(0)),
                comg1 * comg2 + somg1 * somg2);
  k2 = sq(calp0) * g->ep2;
  eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
  C3f(g, eps, C3a);
  B312 = (SinCosSeries(TRUE, ssig2, csig2, C3a, nC3-1) -
          SinCosSeries(TRUE, ssig1, csig1, C3a, nC3-1));
  h0 = -g->f * A3f(g, eps);
  domg12 = salp0 * h0 * (sig12 + B312);
  lam12 = omg12 + domg12;

  if (diffp) {
    if (calp2 == 0)
      dlam12 = - 2 * g->f1 * dn1 / sbet1;
    else {
      real dummy;
      Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
              cbet1, cbet2, &dummy, &dlam12, &dummy,
              FALSE, &dummy, &dummy, C1a, C2a);
      dlam12 *= g->f1 / (calp2 * cbet2);
    }
  }

  *psalp2 = salp2;
  *pcalp2 = calp2;
  *psig12 = sig12;
  *pssig1 = ssig1;
  *pcsig1 = csig1;
  *pssig2 = ssig2;
  *pcsig2 = csig2;
  *peps = eps;
  *pdomg12 = domg12;
  if (diffp)
    *pdlam12 = dlam12;

  return lam12;
}

real A3f(const struct Geodesic* g, real eps) {
  /* Evaluate sum(A3x[k] * eps^k, k, 0, nA3x-1) by Horner's method */
  real v = 0;
  int i;
  for (i = nA3x; i; )
    v = eps * v + g->A3x[--i];
  return v;
}

void C3f(const struct Geodesic* g, real eps, real c[]) {
  /* Evaluate C3 coeffs by Horner's method
   * Elements c[1] thru c[nC3 - 1] are set */
  int i, j, k;
  real mult = 1;
  for (j = nC3x, k = nC3 - 1; k; ) {
    real t = 0;
    for (i = nC3 - k; i; --i)
      t = eps * t + g->C3x[--j];
    c[k--] = t;
  }

  for (k = 1; k < nC3; ) {
    mult *= eps;
    c[k++] *= mult;
  }
}

void C4f(const struct Geodesic* g, real eps, real c[]) {
  /* Evaluate C4 coeffs by Horner's method
   * Elements c[0] thru c[nC4 - 1] are set */
  int i, j, k;
  real mult = 1;
  for (j = nC4x, k = nC4; k; ) {
    real t = 0;
    for (i = nC4 - k + 1; i; --i)
      t = eps * t + g->C4x[--j];
    c[--k] = t;
  }

  for (k = 1; k < nC4; ) {
    mult *= eps;
    c[k++] *= mult;
  }
}

/* Generated by Maxima on 2010-09-04 10:26:17-04:00 */

/* The scale factor A1-1 = mean value of (d/dsigma)I1 - 1 */
real A1m1f(real eps)  {
  real
    eps2 = sq(eps),
    t = eps2*(eps2*(eps2+4)+64)/256;
  return (t + eps) / (1 - eps);
}

/* The coefficients C1[l] in the Fourier expansion of B1 */
void C1f(real eps, real c[])  {
  real
    eps2 = sq(eps),
    d = eps;
  c[1] = d*((6-eps2)*eps2-16)/32;
  d *= eps;
  c[2] = d*((64-9*eps2)*eps2-128)/2048;
  d *= eps;
  c[3] = d*(9*eps2-16)/768;
  d *= eps;
  c[4] = d*(3*eps2-5)/512;
  d *= eps;
  c[5] = -7*d/1280;
  d *= eps;
  c[6] = -7*d/2048;
}

/* The coefficients C1p[l] in the Fourier expansion of B1p */
void C1pf(real eps, real c[])  {
  real
    eps2 = sq(eps),
    d = eps;
  c[1] = d*(eps2*(205*eps2-432)+768)/1536;
  d *= eps;
  c[2] = d*(eps2*(4005*eps2-4736)+3840)/12288;
  d *= eps;
  c[3] = d*(116-225*eps2)/384;
  d *= eps;
  c[4] = d*(2695-7173*eps2)/7680;
  d *= eps;
  c[5] = 3467*d/7680;
  d *= eps;
  c[6] = 38081*d/61440;
}

/* The scale factor A2-1 = mean value of (d/dsigma)I2 - 1 */
real A2m1f(real eps)  {
  real
    eps2 = sq(eps),
    t = eps2*(eps2*(25*eps2+36)+64)/256;
  return t * (1 - eps) - eps;
}

/* The coefficients C2[l] in the Fourier expansion of B2 */
void C2f(real eps, real c[])  {
  real
    eps2 = sq(eps),
    d = eps;
  c[1] = d*(eps2*(eps2+2)+16)/32;
  d *= eps;
  c[2] = d*(eps2*(35*eps2+64)+384)/2048;
  d *= eps;
  c[3] = d*(15*eps2+80)/768;
  d *= eps;
  c[4] = d*(7*eps2+35)/512;
  d *= eps;
  c[5] = 63*d/1280;
  d *= eps;
  c[6] = 77*d/2048;
}

/* The scale factor A3 = mean value of (d/dsigma)I3 */
void A3coeff(struct Geodesic* g) {
  g->A3x[0] = 1;
  g->A3x[1] = (g->n-1)/2;
  g->A3x[2] = (g->n*(3*g->n-1)-2)/8;
  g->A3x[3] = ((-g->n-3)*g->n-1)/16;
  g->A3x[4] = (-2*g->n-3)/64;
  g->A3x[5] = -3/(real)(128);
}

/* The coefficients C3[l] in the Fourier expansion of B3 */
void C3coeff(struct Geodesic* g) {
  g->C3x[0] = (1-g->n)/4;
  g->C3x[1] = (1-g->n*g->n)/8;
  g->C3x[2] = ((3-g->n)*g->n+3)/64;
  g->C3x[3] = (2*g->n+5)/128;
  g->C3x[4] = 3/(real)(128);
  g->C3x[5] = ((g->n-3)*g->n+2)/32;
  g->C3x[6] = ((-3*g->n-2)*g->n+3)/64;
  g->C3x[7] = (g->n+3)/128;
  g->C3x[8] = 5/(real)(256);
  g->C3x[9] = (g->n*(5*g->n-9)+5)/192;
  g->C3x[10] = (9-10*g->n)/384;
  g->C3x[11] = 7/(real)(512);
  g->C3x[12] = (7-14*g->n)/512;
  g->C3x[13] = 7/(real)(512);
  g->C3x[14] = 21/(real)(2560);
}

/* Generated by Maxima on 2012-10-19 08:02:34-04:00 */

/* The coefficients C4[l] in the Fourier expansion of I4 */
void C4coeff(struct Geodesic* g) {
  g->C4x[0] = (g->n*(g->n*(g->n*(g->n*(100*g->n+208)+572)+3432)-12012)+30030)/
    45045;
  g->C4x[1] = (g->n*(g->n*(g->n*(64*g->n+624)-4576)+6864)-3003)/15015;
  g->C4x[2] = (g->n*((14144-10656*g->n)*g->n-4576)-858)/45045;
  g->C4x[3] = ((-224*g->n-4784)*g->n+1573)/45045;
  g->C4x[4] = (1088*g->n+156)/45045;
  g->C4x[5] = 97/(real)(15015);
  g->C4x[6] = (g->n*(g->n*((-64*g->n-624)*g->n+4576)-6864)+3003)/135135;
  g->C4x[7] = (g->n*(g->n*(5952*g->n-11648)+9152)-2574)/135135;
  g->C4x[8] = (g->n*(5792*g->n+1040)-1287)/135135;
  g->C4x[9] = (468-2944*g->n)/135135;
  g->C4x[10] = 1/(real)(9009);
  g->C4x[11] = (g->n*((4160-1440*g->n)*g->n-4576)+1716)/225225;
  g->C4x[12] = ((4992-8448*g->n)*g->n-1144)/225225;
  g->C4x[13] = (1856*g->n-936)/225225;
  g->C4x[14] = 8/(real)(10725);
  g->C4x[15] = (g->n*(3584*g->n-3328)+1144)/315315;
  g->C4x[16] = (1024*g->n-208)/105105;
  g->C4x[17] = -136/(real)(63063);
  g->C4x[18] = (832-2560*g->n)/405405;
  g->C4x[19] = -128/(real)(135135);
  g->C4x[20] = 128/(real)(99099);
}
