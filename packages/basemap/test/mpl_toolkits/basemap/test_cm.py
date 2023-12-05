"""Import test for :mod:`mpl_toolkits.basemap.cm`."""

try:
    import unittest2 as unittest
except ImportError:
    import unittest

from mpl_toolkits.basemap import cm


class TestMplToolkitsBasemapCm(unittest.TestCase):
    """Unittest class for :mod:`mpl_toolkits.basemap.cm`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

    def test_cm_contents(self):
        """Test :mod:`mpl_toolkits.basemap.cm` contents."""

        cmaps = ["GMT_drywet", "GMT_gebco", "GMT_globe", "GMT_haxby",
                 "GMT_no_green", "GMT_ocean", "GMT_polar", "GMT_red2green",
                 "GMT_relief", "GMT_seis", "GMT_split", "GMT_wysiwyg",
                 "s3pcpn", "s3pcpn_l", "StepSeq", "sstanom"]

        self.assertEqual(len(cm.datad), 2 * len(cmaps))
        for cmap in cmaps:
            self.assertTrue(hasattr(cm, cmap))
            self.assertTrue(hasattr(cm, "{0}_r".format(cmap)))


if __name__ == "__main__":
    unittest.main()
