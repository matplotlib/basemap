from __future__ import (absolute_import, division, print_function)

from distutils.version import LooseVersion
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import pyproj

# beginnings of a test suite.

from numpy.testing import TestCase, assert_almost_equal

try:
	from unittest import skipIf
except ImportError:
	# for new features, fallback to unittest backport for Python 2.4 - 2.6
	from unittest2 import skipIf

# For Python 3.x this will be true
PY3 = (sys.version_info[0] == 3)

class TestRotateVector(TestCase):

    def make_array(self):
        lat = np.array([0, 45, 75, 90])
        lon = np.array([0,90,180,270])
        u = np.ones((len(lat), len(lon)))
        v = np.zeros((len(lat), len(lon)))
        return u,v,lat,lon

    def test_cylindrical(self):
        # Cylindrical case
        B = Basemap()
        u,v,lat,lon=self.make_array()
        ru, rv = B.rotate_vector(u,v, lon, lat)
        # Check that the vectors are identical.
        assert_almost_equal(ru, u)
        assert_almost_equal(rv, v)

    def test_nan(self):
        B = Basemap()
        u,v,lat,lon=self.make_array()
        # Set one element to 0, so that the vector magnitude is 0.
        u[1,1] = 0.
        ru, rv = B.rotate_vector(u,v, lon, lat)
        assert not np.isnan(ru).any()
        assert_almost_equal(u, ru)
        assert_almost_equal(v, rv)

    def test_npstere(self):
        # NP Stereographic case
        B=Basemap(projection='npstere', boundinglat=50., lon_0=0.)
        u,v,lat,lon=self.make_array()
        v = np.ones((len(lat), len(lon)))
        ru, rv = B.rotate_vector(u,v, lon, lat)
        assert_almost_equal(ru[2, :],[1,-1,-1,1], 6)
        assert_almost_equal(rv[2, :],[1,1,-1,-1], 6)

class TestShiftGrid(TestCase):

    def make_data_cyc(self):
        loncyc  =  np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300,\
            330, 360],dtype=np.float)
        gridcyc = np.array([[0,  1,  2,  3,   4,   5,   6,   7,   8,   9,  10,\
            11,   0]],dtype=np.float)
        lonoutcyc  =  np.array([-180, -150, -120, -90, -60, -30, 0, 30,60,90,\
            120, 150, 180],dtype=np.float)
        gridoutcyc = np.array([[   6,    7,   8,    9,  10,  11, 0,  1,  2,3,\
            4,   5,   6]],dtype=np.float)
        return loncyc, gridcyc, lonoutcyc, gridoutcyc

    def make_data_nocyc(self):
        lonnocyc  =  np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270,\
            300, 330],dtype=np.float)
        gridnocyc = np.array([[0,  1,  2,  3,   4,   5,   6,   7,   8,   9,\
            10,  11]],dtype=np.float)
        lonoutnocyc  =  np.array([-180, -150, -120, -90, -60, -30, 0, 30, 60,\
            90, 120, 150],dtype=np.float)
        gridoutnocyc = np.array([[   6,    7,   8,    9,  10,  11, 0,  1,  2,\
            3,   4,   5]],dtype=np.float)
        return lonnocyc, gridnocyc, lonoutnocyc, gridoutnocyc

    def make_data_nocyc2(self):
        lonnocyc2  =  np.array([15, 45, 75, 105, 135, 165, 195, 225, 255, 285,\
            315, 345],dtype=np.float)
        gridnocyc2 = np.array([[0,  1,  2,  3,   4,   5,   6,   7,   8,   9,\
            10,  11]],dtype=np.float)
        lonoutnocyc2  =  np.array([-165, -135, -105, -75, -45, -15, 15,45,75,\
            105, 135, 165],dtype=np.float)
        gridoutnocyc2 = np.array([[   6,    7,   8,    9,  10,  11, 0,  1,  2,\
            3,   4,   5]],dtype=np.float)
        return lonnocyc2, gridnocyc2, lonoutnocyc2, gridoutnocyc2

    def test_cyc(self):
        lonin, gridin, lonout, gridout = self.make_data_cyc()
        grid, lon = shiftgrid(lonin[len(lonin)//2], gridin, lonin, start=False)
        assert (lon==lonout).all()
        assert (grid==gridout).all()

    def test_no_cyc(self):
        lonin, gridin, lonout, gridout = self.make_data_nocyc()
        grid, lon = shiftgrid(lonin[len(lonin)//2], gridin, lonin, start=False)
        assert (lon==lonout).all()
        assert (grid==gridout).all()

    def test_no_cyc2(self):
        lonin, gridin, lonout, gridout = self.make_data_nocyc2()
        grid, lon = shiftgrid(lonin[len(lonin)//2], gridin, lonin, start=False)
        assert (lon==lonout).all()
        assert (grid==gridout).all()


class TestShiftdata(TestCase):

    def _get_2d_lons(self, lons1d):
        """
        Generate a 2d grid
        """
        lats = [10, ] * len(lons1d)
        return np.meshgrid(lons1d, lats)[0]

    def test_non_monotonous_longitudes(self):
        """
        when called for scatter, the longitudes passed to shiftdata are
        not necessarily monotonous...
        """
        lons = [179, 180, 180, 0, 290, 10, 320, -150, 350, -250, 250]
        bm = Basemap(lon_0=0)

        # before, having several break points would cause the exception,
        # inside the shiftdata method called from scatter method.
        self.assertRaises(ValueError, bm.shiftdata, lons, fix_wrap_around=True)

        lons_new = bm.shiftdata(lons, fix_wrap_around=False)

        # Check if the modified longitudes are inside of the projection region
        for lon in lons_new:
            assert lon >= bm.projparams["lon_0"] - 180
            assert lon <= bm.projparams["lon_0"] + 180


    def test_shiftdata_on_monotonous_lons(self):
        """
        Test that shiftdata with fix_wrap_around keyword added works as before,
        when it is True
        """

        bm = Basemap(lon_0=0)

        lons_in = [120, 140, 160, 180, 200, 220]
        lons_out_expect = [-160, -140,  120,  140,  160,  180]
        lons_out = bm.shiftdata(lons_in, fix_wrap_around=True)

        assert_almost_equal(lons_out, lons_out_expect)



    def test_2_points_should_work(self):
        """
        Shiftdata should work with 2 points
        """
        bm = Basemap(llcrnrlon=0, llcrnrlat=-80, urcrnrlon=360, urcrnrlat=80, projection='mill')

        lons_expected = [10, 15, 20]
        lonsout = bm.shiftdata(lons_expected[:])
        assert_almost_equal(lons_expected, lonsout)

        lonsout_expected = bm.shiftdata([10, 361, 362])
        lonsout = bm.shiftdata([10, 361])
        assert_almost_equal(lonsout_expected[:len(lonsout)], lonsout)

    def test_1_point_should_work(self):
        bm = Basemap(llcrnrlon=0, llcrnrlat=-80, urcrnrlon=360, urcrnrlat=80, projection='mill')

        # should not fail
        lonsout = bm.shiftdata([361])
        assert_almost_equal(lonsout, [1.0,])

        lonsout = bm.shiftdata([10])
        assert_almost_equal(lonsout, [10.0,])

        lonsin = np.array([361.0])
        lonsin.shape = (1, 1)
        lonsout = bm.shiftdata(lonsin[:])
        assert_almost_equal(lonsout.squeeze(), [1.0,])

    def test_less_than_n_by_3_points_should_work(self):
        bm = Basemap(llcrnrlon=0, llcrnrlat=-80, urcrnrlon=360, urcrnrlat=80, projection='mill')
        lons_expected = self._get_2d_lons([10, 15, 20])

        # nothing should change
        lonsout = bm.shiftdata(lons_expected)
        assert_almost_equal(lons_expected, lonsout)

        # shift n x 3 and n x 2 grids and compare results over overlapping region
        lonsin = self._get_2d_lons([10, 361, 362])
        lonsout_expected = bm.shiftdata(lonsin[:])[:, :2]
        lonsout = bm.shiftdata(lonsin[:, :2])
        assert_almost_equal(lonsout_expected, lonsout)

@skipIf(PY3 and LooseVersion(pyproj.__version__) <= LooseVersion("1.9.4"),
        "Test skipped in Python 3.x with pyproj version 1.9.4 and below.")
class TestProjectCoords(TestCase):
    def get_data(self):
        lons, lats = np.arange(-180, 180, 20), np.arange(-90, 90, 10)
        lats, lons = np.meshgrid(lats, lons)
        lons, lats = lons.copy(order="F"), lats.copy(order="F")
        return lons, lats, Basemap(projection="sinu", lon_0=0)


    def test_convert(self):
        """
        Should not fail on C non-contiguous arrays
        """
        lons, lats, bmp = self.get_data()
        assert not lons.flags['C_CONTIGUOUS']
        assert isinstance(lons, np.ndarray)
        assert isinstance(bmp, Basemap)

        xx1, yy1 = bmp(lons, lats)


    def test_results_should_be_same_for_c_and_f_order_arrays(self):
        lons, lats, bmp = self.get_data()

        xx1, yy1 = bmp(lons.copy(order="C"), lats.copy(order="C"))
        xx2, yy2 = bmp(lons.copy(order="F"), lats.copy(order="F"))

        assert_almost_equal(xx1, xx2)
        assert_almost_equal(yy1, yy2)


class TestInputValidation(TestCase):
    def test_optional_casting(self):
        # Test for the bug reported in gh:#260
        d = {'llcrnrlat': 28.979408, 'urcrnrlat': 35.19622,
             'llcrnrlon': -95.614105, 'urcrnrlon': -77.554749,
             'lon_0': -87.0, 'resolution': 'c', 'lat_0': 32.070374,
             'projection': 'lcc'}
        bmap1 = Basemap(lat_1=30.0, **d)
        bmap2 = Basemap(lat_1=np.array([30.0], dtype='float32'), **d)
        assert bmap1.proj4string == bmap2.proj4string


class TestOrthoProjPolygons(TestCase):
    def test_basemapcreation_should_not_fail(self):
        # different resolutions should work
        for r in ['c', 'l', 'i', 'h', 'f']:
            m = Basemap(projection='ortho',resolution=r,lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)
        pass

@skipIf(not PY3, "Test skipped for Python 2.x as it requires PIL installed")
class TestArcgisimage(TestCase):
	def test_cyl_proj_should_not_fail(self):
		m = Basemap(projection='cyl',
					llcrnrlon=-90,llcrnrlat=30,
					urcrnrlon=-60,urcrnrlat=60)
		m.arcgisimage(verbose=True)



def test():
    """
    Run some tests.
    """
    import unittest
    from . import test
    runner = unittest.TextTestRunner()
    suite = unittest.findTestCases(test)
    runner.run(suite)


if __name__ == '__main__':
    # When called with the -v / --verbose commandline parameter, it will
    # give package dependent version information:
    #   $ python test.py --verbose

    import sys
    import unittest

    from mpl_toolkits.basemap.diagnostic import package_versions

    if '--verbose' in sys.argv or '-v' in sys.argv:
        pkg_vers = package_versions()
        print('Basemaps installed package versions:')
        print('{0}\n'.format(pkg_vers))

    unittest.main()
