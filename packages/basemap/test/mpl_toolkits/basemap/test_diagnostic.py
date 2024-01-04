"""Import test for :mod:`mpl_toolkits.basemap.diagnostic`."""

from collections import namedtuple
try:
    import unittest2 as unittest
except ImportError:
    import unittest

from mpl_toolkits.basemap import diagnostic


class TestMplToolkitsBasemapDiagnostic(unittest.TestCase):
    """Unittest class for :mod:`mpl_toolkits.basemap.diagnostic`."""

    def setUp(self):
        """Define the setup of test scope variables."""

    def tearDown(self):
        """Define the teardown of test scope variables."""

    def test_proj4_version(self):
        """Test getting PROJ version through :mod:`pyproj`."""

        proj_version = diagnostic.proj4_version()
        self.assertIsInstance(proj_version, str)

    def test_package_versions(self):
        """Test getting versions for package dependencies."""

        dependencies = diagnostic.package_versions()
        self.assertIsInstance(dependencies, tuple)
        self.assertEqual(len(dependencies), 11)

    def test_check_proj_inv_hammer(self):
        """Test check for inverse of Hammer project support by PROJ."""

        result = diagnostic.check_proj_inv_hammer(segfault_protection=True)
        self.assertIn(result, [True, False, "Unknown"])


if __name__ == "__main__":
    unittest.main()
