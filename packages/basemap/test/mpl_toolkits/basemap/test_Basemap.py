"""Import test for the :mod:`mpl_toolkits.basemap.Basemap` class."""

import os
import shutil
import tempfile
import datetime as dt
try:
    import unittest2 as unittest
except ImportError:
    import unittest

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.contour import QuadContourSet
from matplotlib.image import AxesImage
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid

try:
    import PIL
except ImportError:
    PIL = None


mpl_version = tuple(map(int, mpl.__version__.split(".")[:2]))


class TestMplToolkitsBasemapBasemap(unittest.TestCase):
    """Unittest class for the :mod:`mpl_toolkits.basemap.Basemap` class."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

        axs_obj = plt.gca()
        axs_obj.clear()

    def test_init_with_ortho(self, resolution=None):
        """Test init with orthographic projection and a given resolution."""

        bmap = Basemap(projection="ortho", resolution=resolution,
                       lat_1=45, lat_2=55, lat_0=50, lon_0=-107)
        self.assertIsInstance(bmap, Basemap)

    def test_init_with_ortho_c_resolution(self):
        """Test init with orthographic projection and a given resolution."""

        self.test_init_with_ortho(resolution="c")

    def test_init_with_ortho_l_resolution(self):
        """Test init with orthographic projection and a given resolution."""

        self.test_init_with_ortho(resolution="l")

    def test_init_with_ortho_i_resolution(self):
        """Test init with orthographic projection and a given resolution."""

        self.test_init_with_ortho(resolution="i")

    def test_init_with_optional_casting(self):
        """Test init for the bug reported in GitHub issue #260."""

        kwds = {
            "llcrnrlat":
                28.979408,
            "urcrnrlat":
                35.19622,
            "llcrnrlon":
                -95.614105,
            "urcrnrlon":
                -77.554749,
            "lon_0":
                -87.0,
            "resolution":
                "c",
            "lat_0":
                32.070374,
            "projection":
                "lcc"
        }

        bmap1_lat_1 = 30.0
        bmap2_lat_1 = np.array([30.0], dtype=np.float32)

        bmap1 = Basemap(lat_1=bmap1_lat_1, **kwds)
        bmap2 = Basemap(lat_1=bmap2_lat_1, **kwds)
        self.assertEqual(bmap1.proj4string, bmap2.proj4string)

    def test_drawcoastlines(self, axs=None, axslen0=10):
        """Test that no lines are missing when drawing coastlines."""

        axs_obj = plt.gca() if axs is None else axs
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0)

        bmap = Basemap(projection="merc", resolution="i", lat_ts=20,
                       llcrnrlat=36.0, llcrnrlon=6.0,
                       urcrnrlat=47.7, urcrnrlon=19.0)

        collection = bmap.drawcoastlines(linewidth=1, color="red")
        self.assertIsInstance(collection, LineCollection)

        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0 + 1)

        lines = collection.get_paths()
        self.assertEqual(len(lines), 27)

    def test_drawcountries(self, axs=None, axslen0=10):
        """Test that no lines are missing when drawing country boundaries."""

        axs_obj = plt.gca() if axs is None else axs
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0)

        bmap = Basemap(projection="merc", resolution="i", lat_ts=20,
                       llcrnrlat=36.0, llcrnrlon=6.0,
                       urcrnrlat=47.7, urcrnrlon=19.0)

        collection = bmap.drawcountries(linewidth=1, color="blue")
        self.assertIsInstance(collection, LineCollection)

        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0 + 1)

        lines = collection.get_paths()
        self.assertEqual(len(lines), 29)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_arcgisimage_with_cyl(self, axs=None, axslen0=10):
        """Test showing an ArcGIS image as background."""

        axs_obj = plt.gca() if axs is None else axs
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0)

        bmap = Basemap(ax=axs, projection="cyl", resolution=None,
                       llcrnrlon=-90, llcrnrlat=30,
                       urcrnrlon=-60, urcrnrlat=60)
        img = bmap.arcgisimage(verbose=False)
        self.assertIsInstance(img, AxesImage)

        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0 + 1)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_arcgisimage_with_cyl_using_cache(self, existing=False, axs=None, axslen0=10):
        """Test showing an ArcGIS image as background."""

        axs_obj = plt.gca() if axs is None else axs
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0)

        bmap = Basemap(ax=axs, projection="cyl", resolution=None,
                       llcrnrlon=-90, llcrnrlat=30,
                       urcrnrlon=-60, urcrnrlat=60)

        # Create cache directory string and check it is empty.
        tmpdir = tempfile.mkdtemp(prefix="tmp-basemap-cachedir-")
        cachedir = tmpdir if existing else os.path.join(tmpdir, "cachedir")
        if os.path.isdir(cachedir):
            self.assertEqual(len(os.listdir(cachedir)), 0)

        try:
            # Check that the first call populates the cache.
            img = bmap.arcgisimage(verbose=False, cachedir=cachedir)
            self.assertEqual(len(os.listdir(cachedir)), 1)
            # Check output properties after the first call.
            self.assertIsInstance(img, AxesImage)
            axs_children = axs_obj.get_children()
            self.assertEqual(len(axs_children), axslen0 + 1)
            # Check that the second call does not update the cache.
            img = bmap.arcgisimage(verbose=False, cachedir=cachedir)
            self.assertEqual(len(os.listdir(cachedir)), 1)
            # Check output properties after the second call.
            self.assertIsInstance(img, AxesImage)
            axs_children = axs_obj.get_children()
            self.assertEqual(len(axs_children), axslen0 + 2)
        finally:
            if os.path.isdir(tmpdir):
                shutil.rmtree(tmpdir)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_arcgisimage_with_cyl_using_cache_already_existing(self):
        """Test showing an ArcGIS image as background."""

        self.test_arcgisimage_with_cyl_using_cache(existing=True)

    def _test_basemap_data_warpimage(self, method, axs=None, axslen0=10):
        """Test drawing a map background from :mod:`basemap_data`."""

        axs_obj = plt.gca() if axs is None else axs
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0)

        bmap = Basemap(ax=axs, projection="moll", resolution=None, lon_0=0)
        img = getattr(bmap, method)(ax=axs, scale=0.1)
        self.assertIsInstance(img, AxesImage)

        flag = int(mpl_version < (3, 5))
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0 + 3)
        self.assertIsInstance(axs_children[1 - flag], Polygon)
        self.assertIsInstance(axs_children[2 - flag], Polygon)
        self.assertIsInstance(axs_children[(axslen0 + 1) * flag], AxesImage)
        self.assertIs(axs_children[(axslen0 + 1) * flag], img)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_bluemarble(self, axs=None, axslen0=10):
        """Test drawing a map with a blue marble image as background."""

        method = "bluemarble"
        self._test_basemap_data_warpimage(method, axs=axs, axslen0=axslen0)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_bluemarble_with_custom_axes(self):
        """Test drawing a map with a blue marble image as background."""

        _, axs = plt.subplots()
        self.test_bluemarble(axs=axs, axslen0=10)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_etopo(self, axs=None, axslen0=10):
        """Test drawing a map with an ETOPO relief image as background."""

        method = "etopo"
        self._test_basemap_data_warpimage(method, axs=axs, axslen0=axslen0)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_etopo_with_custom_axes(self):
        """Test drawing a map with an ETOPO relief image as background."""

        _, axs = plt.subplots()
        self.test_etopo(axs=axs, axslen0=10)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_shadedrelief(self, axs=None, axslen0=10):
        """Test drawing a map with a shaded relief image as background."""

        method = "shadedrelief"
        self._test_basemap_data_warpimage(method, axs=axs, axslen0=axslen0)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_shadedrelief_with_custom_axes(self):
        """Test drawing a map with a shaded relief image as background."""

        _, axs = plt.subplots()
        self.test_shadedrelief(axs=axs, axslen0=10)

    def _test_generic_contour_function(self, function):
        """Generic test for the `contour` and `contourf` methods."""

        bmap = Basemap(projection="ortho", lat_0=45, lon_0=-100, resolution=None)

        # Create a regular lat/lon grid.
        nlats = 73
        nlons = 145
        delta = 2 * np.pi / (nlons - 1)
        indx = np.indices((nlats, nlons))
        lats = (0.5 * np.pi - delta * indx[0, :, :])
        lons = (delta * indx[1, :, :])

        # Create some data the regular lat/lon grid.
        mean = 0.50 * np.cos(2 * lats) * ((np.sin(2 * lats))**2 + 2)
        wave = 0.75 * np.cos(4 * lons) * np.sin(2 * lats)**8
        data = mean + wave

        # Compute native map projection coordinates of lat/lon grid.
        x, y = bmap(np.degrees(lons), np.degrees(lats))

        # Contour data over the map and check output.
        cset = getattr(bmap, function)(x, y, data, 15)
        self.assertIsInstance(cset, QuadContourSet)

    def test_contour(self):
        """Test drawing contours on a map."""

        self._test_generic_contour_function("contour")

    def test_contourf(self):
        """Test drawing filled contours on a map."""

        self._test_generic_contour_function("contourf")

    def test_nightshade(self):
        """Test drawing the day/night terminator and night shade on a map."""

        bmap = Basemap(projection="mill", lon_0=180)
        cset = bmap.nightshade(date=dt.datetime(1970, 1, 1))
        self.assertIsInstance(cset, QuadContourSet)


class TestMplToolkitsBasemapBasemapCall(unittest.TestCase):
    """Unittest class for :meth:`mpl_toolkits.basemap.Basemap.__call__`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

        axs_obj = plt.gca()
        axs_obj.clear()

    @staticmethod
    def get_input_objects():
        """Return geographic coordinate arrays and :class:`Basemap` instance."""

        lons, lats = np.arange(-180, 180, 20), np.arange(-90, 90, 10)
        lats, lons = np.meshgrid(lats, lons)
        lons, lats = lons.copy(order="F"), lats.copy(order="F")
        return Basemap(projection="sinu", lon_0=0), lons, lats

    def test_transform_with_c_non_contiguous_arrays(self):
        """Test transform with C non-contiguous arrays."""

        bmap, lons, lats = self.get_input_objects()
        self.assertIsInstance(bmap, Basemap)
        self.assertIsInstance(lats, np.ndarray)
        self.assertIsInstance(lons, np.ndarray)
        self.assertFalse(lats.flags["C_CONTIGUOUS"])
        self.assertFalse(lons.flags["C_CONTIGUOUS"])

        xx1, yy1 = bmap(lons, lats)
        self.assertIsInstance(xx1, np.ndarray)
        self.assertIsInstance(yy1, np.ndarray)
        self.assertEqual(xx1.shape, lons.shape)
        self.assertEqual(yy1.shape, lats.shape)

    def test_transform_equal_with_c_and_f_order_arrays(self):
        """Test transform with C contiguous arrays and F contiguous arrays."""

        bmap, lons, lats = self.get_input_objects()

        xx1, yy1 = bmap(lons.copy(order="C"), lats.copy(order="C"))
        xx2, yy2 = bmap(lons.copy(order="F"), lats.copy(order="F"))
        self.assertTrue(np.allclose(xx1, xx2))
        self.assertTrue(np.allclose(yy1, yy2))


class TestMplToolkitsBasemapBasemapShiftData(unittest.TestCase):
    """Unittest class for :meth:`mpl_toolkits.basemap.Basemap.shiftdata`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

        axs_obj = plt.gca()
        axs_obj.clear()

    def generate_2d_longitude_grid(self, lons1d):
        """Return a 2D longitude grid."""

        lats = [10] * len(lons1d)
        return np.meshgrid(lons1d, lats)[0]

    def test_shiftdata_with_non_monotonous_longitudes(self):
        """Test that shiftdata works with non-monontonous longitudes.

        For example, this is the case with scatter data.
        """

        # Before, having several break points would cause the exception
        # inside the `shiftdata` method called from `scatter` method.
        lons = [179, 180, 180, 0, 290, 10, 320, -150, 350, -250, 250]
        bmap = Basemap(projection="cyl", resolution=None, lon_0=0)
        lons_new = bmap.shiftdata(lons, fix_wrap_around=True)
        for lon in lons_new:
            self.assertGreaterEqual(lon, bmap.projparams["lon_0"] - 180)
            self.assertLessEqual(lon, bmap.projparams["lon_0"] + 180)

        # Check if the modified longitudes are inside the projection region.
        lons_new = bmap.shiftdata(lons, fix_wrap_around=False)
        for lon in lons_new:
            self.assertGreaterEqual(lon, bmap.projparams["lon_0"] - 180)
            self.assertLessEqual(lon, bmap.projparams["lon_0"] + 180)

    def test_shiftdata_with_monotonous_lons(self):
        """Test that shiftdata works with `fix_wrap_around=True` as before."""

        bmap = Basemap(projection="cyl", resolution=None, lon_0=0)

        lons_in = [120, 140, 160, 180, 200, 220]
        lons_out_expected = [-160, -140, 120, 140, 160, 180]

        lons_out = bmap.shiftdata(lons_in, fix_wrap_around=True)
        self.assertTrue(np.allclose(lons_out, lons_out_expected))

    def test_shiftdata_with_1_point(self):
        """Test that shiftdata works with 1 points."""

        bmap = Basemap(projection="mill", resolution=None,
                       llcrnrlon=0, llcrnrlat=-80,
                       urcrnrlon=360, urcrnrlat=80)

        # This should not fail due to longitude out of bounds.
        lons_out = bmap.shiftdata([361])
        self.assertTrue(np.allclose(lons_out, [1.0]))

        lons_out = bmap.shiftdata([10])
        self.assertTrue(np.allclose(lons_out, [10.0]))

        lons_in = np.asarray([[361.0]])
        lons_out = bmap.shiftdata(lons_in[:])
        self.assertTrue(np.allclose(lons_out, [[1.0]]))

    def test_shiftdata_with_2_points(self):
        """Test that shiftdata works with 2 points."""

        bmap = Basemap(projection="mill", resolution=None,
                       llcrnrlon=0, llcrnrlat=-80,
                       urcrnrlon=360, urcrnrlat=80)

        lons_out_expected = [10, 15, 20]
        lons_out = bmap.shiftdata(lons_out_expected[:])
        self.assertTrue(np.allclose(lons_out, lons_out_expected))

        lons_out_expected = bmap.shiftdata([10, 361, 362])
        lons_out = bmap.shiftdata([10, 361])
        self.assertTrue(np.allclose(lons_out, lons_out_expected[:lons_out.size]))

    def test_shiftdata_with_less_than_n_by_3_points(self):
        """Test that shiftdata works with `(n x 3)` and `(n x 2)` grids."""

        bmap = Basemap(projection="mill", resolution=None,
                       llcrnrlon=0, llcrnrlat=-80,
                       urcrnrlon=360, urcrnrlat=80)

        # Test that nothing should change here.
        lons_out_expected = self.generate_2d_longitude_grid([10, 15, 20])
        lons_out = bmap.shiftdata(lons_out_expected)
        self.assertTrue(np.allclose(lons_out, lons_out_expected))

        # Shift n x 3 and n x 2 grids and compare results in overlapping area.
        lons_in = self.generate_2d_longitude_grid([10, 361, 362])
        lons_out = bmap.shiftdata(lons_in[:, :2])
        lons_out_expected = bmap.shiftdata(lons_in)[:, :2]
        self.assertTrue(np.allclose(lons_out, lons_out_expected))


class TestMplToolkitsBasemapBasemapRotateVector(unittest.TestCase):
    """Unittest class for :meth:`mpl_toolkits.basemap.Basemap.rotate_vector`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

        axs_obj = plt.gca()
        axs_obj.clear()

    @staticmethod
    def get_input_objects():
        """Return geographic coordinates and vector components."""

        lat = np.array([0, 45, 75, 90])
        lon = np.array([0, 90, 180, 270])

        shape = (lat.size, lon.size)
        u = np.ones(shape)
        v = np.zeros(shape)

        return u, v, lat, lon

    def test_rotate_vector_cylindrical(self):
        """Test vector rotation with cylindrical projection."""

        u, v, lat, lon = self.get_input_objects()

        bmap = Basemap(projection="cyl", resolution=None)

        # Check that the vectors are identical after rotation.
        # pylint: disable=unbalanced-tuple-unpacking
        rotu, rotv = bmap.rotate_vector(u, v, lon, lat, returnxy=False)
        self.assertTrue(np.allclose(rotu, u))
        self.assertTrue(np.allclose(rotv, v))

    def test_rotate_vector_cylindrical_without_nan(self):
        """Test vector rotation with cylindrical projection."""

        # Set one `u` element to 0, so that the vector magnitude is 0.
        u, v, lat, lon = self.get_input_objects()
        u[1, 1] = 0

        bmap = Basemap(projection="cyl", resolution=None)

        # Check that the vectors are identical after rotation.
        # pylint: disable=unbalanced-tuple-unpacking
        rotu, rotv = bmap.rotate_vector(u, v, lon, lat, returnxy=False)
        self.assertFalse(np.isnan(rotu).any())
        self.assertTrue(np.allclose(rotu, u))
        self.assertTrue(np.allclose(rotv, v))

    def test_rotate_vector_npstere(self):
        """Test vector rotation with NP stereographic projection."""

        # Set all `v` elements to 1.
        u, v, lat, lon = self.get_input_objects()
        v = np.ones(v.shape)

        bmap = Basemap(projection="npstere", resolution=None,
                       boundinglat=50., lon_0=0.)

        # pylint: disable=unbalanced-tuple-unpacking
        rotu, rotv = bmap.rotate_vector(u, v, lon, lat, returnxy=False)
        self.assertTrue(np.allclose(rotu[2, :], [+1, -1, -1, +1]))
        self.assertTrue(np.allclose(rotv[2, :], [+1, +1, -1, -1]))


class TestMplToolkitsBasemapShiftGrid(unittest.TestCase):
    """Unittest class for :func:`mpl_toolkits.basemap.shiftgrid`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

    def test_shifgrid_with_cyclic_data(self):
        """Test shiftgrid with some cyclic data."""

        lonin = np.array(
            [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360],
            dtype=np.float64)
        gridin = np.array(
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0]],
            dtype=np.float64)
        lonout = np.array(
            [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180],
            dtype=np.float64)
        gridout = np.array(
            [[6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6]],
            dtype=np.float64)

        grid, lon = shiftgrid(lonin[len(lonin) // 2], gridin, lonin, start=False)
        self.assertTrue((lon == lonout).all())
        self.assertTrue((grid == gridout).all())

    def test_shiftgrid_with_no_cyclicdata_1(self):
        """Test shiftgrid with some no-cyclic data."""

        lonin = np.array(
            [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330],
            dtype=np.float64)
        gridin = np.array(
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]],
            dtype=np.float64)
        lonout = np.array(
            [-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150],
            dtype=np.float64)
        gridout = np.array(
            [[6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5]],
            dtype=np.float64)

        grid, lon = shiftgrid(lonin[len(lonin) // 2], gridin, lonin, start=False)
        self.assertTrue((lon == lonout).all())
        self.assertTrue((grid == gridout).all())

    def test_shiftgrid_with_no_cyclicdata_2(self):
        """Test shiftgrid with some no-cyclic data."""

        lonin = np.array(
            [15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345],
            dtype=np.float64)
        gridin = np.array(
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]],
            dtype=np.float64)
        lonout = np.array(
            [-165, -135, -105, -75, -45, -15, 15, 45, 75, 105, 135, 165],
            dtype=np.float64)
        gridout = np.array(
            [[6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5]],
            dtype=np.float64)

        grid, lon = shiftgrid(lonin[len(lonin) // 2], gridin, lonin, start=False)
        self.assertTrue((lon == lonout).all())
        self.assertTrue((grid == gridout).all())


if __name__ == "__main__":
    unittest.main()
