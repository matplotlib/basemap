"""Import test for :mod:`mpl_toolkits.basemap.proj`."""

try:
    import unittest2 as unittest
except ImportError:
    import unittest

from mpl_toolkits.basemap.proj import Proj


class TestMplToolkitsBasemapProj(unittest.TestCase):
    """Unittest class for :mod:`mpl_toolkits.basemap.proj`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

    def get_awips_example(self):
        """Return :mod:`Proj` instance for AWIPS grid 221 parameters.

        See https://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
        for further information.
        """

        nx, dx = 349, 32463.41
        ny, dy = 277, 32463.41
        projparams = {
            "proj": "lcc",
            "R": 6371200,
            "lat_1": 50,
            "lat_2": 50,
            "lon_0": -107
        }

        llcrnrlon, llcrnrlat = -145.5, 1.0
        urcrnrlon, urcrnrlat = (nx - 1) * dx, (ny - 1) * dy
        awips221 = Proj(projparams, llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
                        urcrnrislatlon=False)

        return awips221

    def test_proj_lcc_xy_lower_left_corner(self):
        """Test `class`:Proj: instace for AWIPS grid 221 parameters."""

        awips221 = self.get_awips_example()

        llcornerx, llcornery = awips221(awips221.llcrnrlon, awips221.llcrnrlat)
        self.assertAlmostEqual(llcornerx, 0)
        self.assertAlmostEqual(llcornery, 0)

    def test_proj_lcc_lonlat_lower_left_corners(self):
        """Test `class`:Proj: instace for AWIPS grid 221 parameters."""

        awips221 = self.get_awips_example()
        llcornerx, llcornery = 0, 0

        llcornerlon, llcornerlat = awips221(llcornerx, llcornery, inverse=True)
        self.assertAlmostEqual(llcornerlon, -145.5, places=3)
        self.assertAlmostEqual(llcornerlat, +1.0, places=3)

    def test_proj_lcc_lonlat_lower_right_corners(self):
        """Test `class`:Proj: instace for AWIPS grid 221 parameters."""

        awips221 = self.get_awips_example()
        lrcornerx, lrcornery = awips221.urcrnrx, 0

        lrcornerlon, lrcornerlat = awips221(lrcornerx, lrcornery, inverse=True)
        self.assertAlmostEqual(lrcornerlon, -68.318, places=3)
        self.assertAlmostEqual(lrcornerlat, +0.897, places=3)

    def test_proj_lcc_lonlat_upper_left_corners(self):
        """Test `class`:Proj: instace for AWIPS grid 221 parameters."""

        awips221 = self.get_awips_example()
        ulcornerx, ulcornery = 0, awips221.urcrnry

        ulcornerlon, ulcornerlat = awips221(ulcornerx, ulcornery, inverse=True)
        self.assertAlmostEqual(ulcornerlon, 148.639, places=3)
        self.assertAlmostEqual(ulcornerlat, 46.635, places=3)

    def test_proj_lcc_lonlat_upper_right_corners(self):
        """Test `class`:Proj: instace for AWIPS grid 221 parameters."""

        awips221 = self.get_awips_example()
        urcornerx, urcornery = awips221.urcrnrx, awips221.urcrnry

        urcornerlon, urcornerlat = awips221(urcornerx, urcornery, inverse=True)
        self.assertAlmostEqual(urcornerlon, -2.566, places=3)
        self.assertAlmostEqual(urcornerlat, 46.352, places=3)

    def test_proj_lcc_makegrid_makegrid3d(self):
        """Test `class`:Proj: instace for AWIPS grid 221 parameters."""

        nx, ny = 349, 277
        awips221 = self.get_awips_example()

        lons, lats = awips221.makegrid(nx, ny, returnxy=False)
        lonlats = awips221.makegrid3d(nx, ny, returnxy=False)

        self.assertTrue((lons == lonlats[..., 0]).all())
        self.assertTrue((lats == lonlats[..., 1]).all())


if __name__ == "__main__":
    unittest.main()
