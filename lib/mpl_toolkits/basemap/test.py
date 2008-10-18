from mpl_toolkits.basemap import Basemap
import numpy as np

# beginnings of a test suite.

from numpy.testing import NumpyTestCase,assert_almost_equal
class TestRotateVector(NumpyTestCase):
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

def test():
    """
    Run some tests.
    """
    import unittest
    suite = unittest.makeSuite(TestRotateVector,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

if __name__ == '__main__':
    test()
