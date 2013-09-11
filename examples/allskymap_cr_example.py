"""
Example of astronomical use of AllSkyMap class in allskymap.py module

Plot an all-sky map showing locations of the 27 highest-energy ultra-high
energy cosmic rays detected by the Auger (south) experiment as of Aug 2007,
and locations of 18 (fictitious!) candidate sources.  Indicate CR direction
uncertainties and source scattering scales with tissots, and show the
nearest candidate source to each CR with geodesics.

Created 2011-02-07 by Tom Loredo
"""

try:
    from cStringIO import StringIO
except:
    from io import StringIO
import numpy as np
from numpy import cos, sin, arccos, deg2rad, rad2deg
import csv, re, sys

import matplotlib.pyplot as plt
from allskymap import AllSkyMap
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
import matplotlib.ticker as ticker


class Source:
    """
    Parse and store data for a celestial source.
    """
    
    int_re = re.compile(r'^[-+]?[0-9]+$')
    # float_re = re.compile(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$')

    def __init__(self, id, year, day, l, b, sig=None, **kwds):
        self.id = int(id)
        self.year = int(year)
        self.day = int(day)
        self.l = float(l)
        self._l = deg2rad(self.l)  # radians
        self.b = float(b)
        self._b = deg2rad(self.b)  # radians
        if sig is not None:
            self.sig = float(sig)
            self._sig = deg2rad(self.sig)
        else:
            self.sig, self._sig = None, None
        # If extra values are specified as keywords, set them as
        # attributes.  The quick way is to use self.__dict__.update(kwds),
        # but we want to convert to int or float values if possible.
        # *** Consider using matplotlib.cbook.converter.
        if kwds:
            for key, val in kwds.items():
                try:
                    nval = float(val)
                    # It's a number; but it may be an int.
                    if self.int_re.match(str(val)):
                        nval = int(val)
                    setattr(self, key, nval)
                except ValueError:  # non-numerical value
                    setattr(self, key, val)

    def gcangle(self, src):
        """
        Calculate the great circle angle to another source.
        """
        # Use the law of cosines; note it is usually expressed with colattitude.
        c = sin(self._b)*sin(src._b) + \
            cos(self._b)*cos(src._b)*cos(self._l - src._l)
        return rad2deg(arccos(c))


# Auger UHE cosmic ray data, Jan 2004 to Aug 2007
# From Appendix A of Abraham et al. (2008); "Correlation of the highest-energy
# cosmic rays with the positions of nearby active galactic nuclei,"
# Astropart.Phys.29:188-204,2008; Erratum-ibid.30:45,2008

# Year	day 	 ang	S(1000)	E (EeV)	RA	Dec	Longitude	Latitude
# * = w/i 3.2 deg of AGN
AugerData = StringIO(
"""2004	125	47.7	252	70	267.1	-11.4	15.4	8.4
2004	142	59.2	212	84	199.7	-34.9	-50.8	27.6	*
2004	282	26.5	328	66	208.0	-60.3	-49.6	1.7	*
2004	339	44.7	316	83	268.5	-61.0	-27.7	-17.0	*
2004	343	23.4	323	63	224.5	-44.2	-34.4	13.0	*
2005	54	35.0	373	84	17.4	-37.9	-75.6	-78.6	*
2005	63	54.5	214	71	331.2	-1.2	58.8	-42.4	*
2005	81	17.2	308	58	199.1	-48.6	-52.8	14.1	*
2005	295	15.4	311	57	332.9	-38.2	4.2	-54.9	*
2005	306	40.1	248	59	315.3	-0.3	48.8	-28.7	*
2005	306	14.2	445	84	114.6	-43.1	-103.7	-10.3
2006	35	30.8	398	85	53.6	-7.8	-165.9	-46.9	*
2006	55	37.9	255	59	267.7	-60.7	-27.6	-16.5	*
2006	81	34.0	357	79	201.1	-55.3	-52.3	7.3
2006	185	59.1	211	83	350.0	9.6	88.8	-47.1	*
2006	296	54.0	208	69	52.8	-4.5	-170.6	-45.7	*
2006	299	26.0	344	69	200.9	-45.3	-51.2	17.2	*
2007	13	14.3	762	148	192.7	-21.0	-57.2	41.8
2007	51	39.2	247	58	331.7	2.9	63.5	-40.2	*
2007	69	30.4	332	70	200.2	-43.4	-51.4	19.2	**
2007	84	17.3	340	64	143.2	-18.3	-109.4	23.8	*
2007	145	23.9	392	78	47.7	-12.8	-163.8	-54.4	*
2007	186	44.8	248	64	219.3	-53.8	-41.7	5.9
2007	193	18.0	469	90	325.5	-33.5	12.1	-49.0	*
2007	221	35.3	318	71	212.7	-3.3	-21.8	54.1	*
2007	234	33.2	365	80	185.4	-27.9	-65.1	34.5
2007	235	42.6	276	69	105.9	-22.9	-125.2	-7.7
""")
AugerTable = csv.reader(AugerData, dialect='excel-tab')
CRs = {}
for id, row in enumerate(AugerTable):
    # Make an integer ID from Year+Day (presumes none on same day!).
    src = Source(id, row[0], row[1], row[7], row[8], E=float(row[4]))
    CRs[src.id] = src
sys.stdout.write('Parsed data for %s UHE CRs...\n'%len(CRs))

# Partly fictitious candidate source locations.
# src.id src.l_deg	src.b_deg	src.xProj	src.yProj
# tab-delimited
CandData = StringIO(
"""1	270.	-28.
2	229.	-80.
3	141.	-47.
4	172.	-51.
5	251.	-51.
6	241.	-36.
7	281.	26.
8	241.	64.
9	240.	64.
10	148.	70.
11	305.	13.
12	98.	79.
13	309.	19.
14	104.	68.
15	104.	68.
16	321.	15.
17	328.	-14.
18	177.5	-35.
""")
# Add this line above to see a tissot overlapping the map limb.
CandTable = csv.reader(CandData, dialect='excel-tab')
cands = {}
for row in CandTable:
    src = Source(row[0], 0, 0, row[1], row[2])
    cands[src.id] = src
sys.stdout.write('Parsed data for %s candidate sources...\n' % len(cands))

# Calculate the separation matrix; track the closest candidate to each CR.
sepn = {}
for cr in CRs.values():
    id, sep = None, 181.
    for cand in cands.values():
        val = cr.gcangle(cand)
        sepn[cr.id, cand.id] = val
        if val < sep:
            id, sep = cand.id, val
    # Store the closest candidate id and separation as a CR attribute.
    cr.near_id = id
    cr.near_ang = sep


# Note that Hammer & Mollweide projections enforce a 2:1 aspect ratio.
# Choose figure size for a 2:1 plot, with room at bottom for colorbar.
fig = plt.figure(figsize=(12,7))
main_ax = plt.axes([0.05, .19, .9, .75])  # rect=L,B,W,H

# Set up the projection and draw a grid.
m = AllSkyMap(ax=main_ax, projection='hammer')
m.drawmapboundary(fill_color='white')
m.drawparallels(np.arange(-75,76,15), linewidth=0.5, dashes=[1,2],
    labels=[1,0,0,0], fontsize=9)
m.drawmeridians(np.arange(-150,151,30), linewidth=0.5, dashes=[1,2])

# Label a subset of meridians.
lons = np.arange(-150,151,30)
m.label_meridians(lons, fontsize=9, vnudge=1,
                halign='left', hnudge=-1)  # nudge<0 shifts to right

# Plot CR directions.
lvals = [src.l for src in CRs.values()]
bvals = [src.b for src in CRs.values()]
x, y = m(lvals, bvals)
# These symbols will be covered by opaque tissots; plot them anyway
# so there is a collection for the legend.
cr_pts = m.scatter(x, y, s=8, c='r', marker='o', linewidths=.5,
    edgecolors='none')

# Plot tissots showing uncertainties, colored by energy.
# We use 1 deg uncertainties, which are probably ~2 sigma for most events.
Evals = np.array([src.E for src in CRs.values()])
norm_E = Normalize(Evals.min()-10, Evals.max()+20)  # -+ for jet_r for brt clrs
# jet_r color map is in spectral sequence, with a small unsaturated
# green range, helping to distinguish CRs from candidates.
cmap = plt.cm.get_cmap('jet_r')
for cr in CRs.values():
    color = cmap(norm_E(cr.E))[0:3]  # ignore alpha
    m.tissot(cr.l, cr.b, 2., 30, ec='none', color=color, alpha=1)

# Plot candidate directions.
lvals = [src.l for src in cands.values()]
bvals = [src.b for src in cands.values()]
x, y = m(lvals, bvals)
cand_pts = m.scatter(x, y, marker='+', linewidths=.5, 
    edgecolors='k', facecolors='none', zorder=10)  # hi zorder -> top

# Plot tissots showing possible scale of candidate scatter.
for l, b in zip(lvals, bvals):
    m.tissot(l, b, 5., 30, ec='none', color='g', alpha=0.25)

# Show the closest candidate to each CR.
for cr in CRs.values():
    cand = cands[cr.near_id]
    m.geodesic(cr.l, cr.b, cand.l, cand.b, lw=0.5, ls='-', c='g')

plt.title('UHE Cosmic Rays and Candidate Sources')
plt.legend([cr_pts, cand_pts], ['UHE CR', 'Candidate'],
    frameon=False, loc='lower right', scatterpoints=1)

# Plot a colorbar for the CR energies.
cb_ax = plt.axes([0.25, .1, .5, .03], frameon=False)  # rect=L,B,W,H
#bar = ColorbarBase(cb_ax, cmap=cmap, orientation='horizontal', drawedges=False)
vals = np.linspace(Evals.min(), Evals.max(), 100)
bar = ColorbarBase(cb_ax, values=vals, norm=norm_E, cmap=cmap, 
    orientation='horizontal', drawedges=False)
bar.set_label('CR Energy (EeV)')

plt.show()
