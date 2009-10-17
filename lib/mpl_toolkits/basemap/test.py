from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np

# beginnings of a test suite.

from numpy.testing import TestCase,assert_almost_equal

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
        grid, lon = shiftgrid(lonin[len(lonin)/2], gridin, lonin, start=False)
        assert (lon==lonout).all()
        assert (grid==gridout).all()

    def test_no_cyc(self):
        lonin, gridin, lonout, gridout = self.make_data_nocyc()
        grid, lon = shiftgrid(lonin[len(lonin)/2], gridin, lonin, start=False)
        assert (lon==lonout).all()
        assert (grid==gridout).all()

    def test_no_cyc2(self):
        lonin, gridin, lonout, gridout = self.make_data_nocyc2()
        grid, lon = shiftgrid(lonin[len(lonin)/2], gridin, lonin, start=False)
        assert (lon==lonout).all()
        assert (grid==gridout).all()


def test():
    """
    Run some tests.
    """
    import unittest
    rotatevector_suite = unittest.makeSuite(TestRotateVector,'test')
    shiftgrid_suite = unittest.makeSuite(TestShiftGrid,'test')
    runner = unittest.TextTestRunner()
    runner.run(rotatevector_suite)
    runner.run(shiftgrid_suite)

if __name__ == '__main__':
    test()
