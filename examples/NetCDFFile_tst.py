import sys
import unittest
import os
import tempfile
from numpy import ma
from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.random.mtrand import uniform 
from mpl_toolkits.basemap import NetCDFFile

# test automatic conversion of masked arrays, and
# packing/unpacking of short ints.

FILE_NAME = tempfile.mktemp(".nc")
ndim = 10
ranarr = 100.*uniform(size=(ndim))
packeddata = 10.*uniform(size=(ndim))
missing_value = -9999.
ranarr[::2] = missing_value
maskedarr = ma.masked_values(ranarr,-9999.)
scale_factor = (packeddata.max()-packeddata.min())/(2.*32766.)
add_offset = 0.5*(packeddata.max()+packeddata.min())
packeddata2 = ((packeddata-add_offset)/scale_factor).astype('i2')

class TestCase(unittest.TestCase):

    def setUp(self):
        self.file = FILE_NAME
        file = NetCDFFile(self.file,'w')
        file.createDimension('n', None) # use unlimited dim.
        foo = file.createVariable('maskeddata', 'f8', ('n',))
        foo.missing_value = missing_value
        bar = file.createVariable('packeddata', 'i2', ('n',))
        bar.scale_factor = scale_factor
        bar.add_offset = add_offset
        foo[0:ndim] = maskedarr
        bar[0:ndim] = packeddata
        file.close()

    def tearDown(self):
        # Remove the temporary files
        os.remove(self.file)

    def runTest(self):
        """testing auto-conversion of masked arrays and packed integers""" 
        # no auto-conversion.
        file = NetCDFFile(self.file,maskandscale=False)
        datamasked = file.variables['maskeddata']
        datapacked = file.variables['packeddata']
        # check missing_value, scale_factor and add_offset attributes.
        assert datamasked.missing_value == missing_value
        assert datapacked.scale_factor == scale_factor
        assert datapacked.add_offset == add_offset
        assert_array_equal(datapacked[:],packeddata2)
        assert_array_almost_equal(datamasked[:],ranarr)
        file.close()
        # auto-conversion
        file = NetCDFFile(self.file)
        datamasked = file.variables['maskeddata']
        datapacked = file.variables['packeddata']
        assert_array_almost_equal(datamasked[:].filled(),ranarr)
        assert_array_almost_equal(datapacked[:],packeddata,decimal=4)
        file.close()

if __name__ == '__main__':
    unittest.main()
