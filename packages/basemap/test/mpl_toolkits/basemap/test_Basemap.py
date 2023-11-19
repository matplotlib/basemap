"""Import test for the :mod:`mpl_toolkits.basemap.Basemap` class."""

try:
    import unittest2 as unittest
except ImportError:
    import unittest

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.image import AxesImage
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

try:
    import PIL
except ImportError:
    PIL = None


mpl_version = tuple(map(int, mpl.__version__.split(".")[:2]))


class TestMplToolkitsBasemapBasemap(unittest.TestCase):
    """Unittest class for the :mod:`mpl_toolkits.basemap.Basemap` class."""

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_shadedrelief(self, axs=None, axslen0=10):
        """Test drawing a map with a shaded relief image as background."""

        axs_obj = plt.gca() if axs is None else axs
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0)

        bmap = Basemap(ax=axs, projection="moll", resolution=None, lon_0=0)
        img = bmap.shadedrelief(ax=axs, scale=0.1)
        self.assertIsInstance(img, AxesImage)

        flag = int(mpl_version < (3, 5))
        axs_children = axs_obj.get_children()
        self.assertEqual(len(axs_children), axslen0 + 3)
        self.assertIsInstance(axs_children[1 - flag], Polygon)
        self.assertIsInstance(axs_children[2 - flag], Polygon)
        self.assertIsInstance(axs_children[(axslen0 + 1) * flag], AxesImage)
        self.assertIs(axs_children[(axslen0 + 1) * flag], img)

    @unittest.skipIf(PIL is None, reason="pillow unavailable")
    def test_shadedrelief_with_custom_axes(self):
        """Test drawing a map with a shaded relief image as background."""

        _, axs = plt.subplots()
        self.test_shadedrelief(axs=axs, axslen0=10)


if __name__ == "__main__":
    unittest.main()
