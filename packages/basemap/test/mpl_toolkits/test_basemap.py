"""Import test for the :mod:`mpl_toolkits.basemap` package."""

try:
    import unittest2 as unittest
except ImportError:
    import unittest

from mpl_toolkits import basemap


class TestMplToolkitsBasemap(unittest.TestCase):
    """Unittest class for the :mod:`mpl_toolkits.basemap` package."""

    def test_version_attribute(self):
        """Test that basic library import is working."""

        self.assertTrue(hasattr(basemap, "__version__"))
        basemap_version = basemap.__version__

        num = r"(?:0|[1-9]\d*)"
        build = r"(?:dev|a[1-4]|b[1-3]|rc[1-2])"
        semver = r"^({0}\.{0}\.{0})(?:[+-]?({1}))?$".format(num, build)
        self.assertRegex(basemap_version, semver)


if __name__ == "__main__":
    unittest.main()
