import numpy as np
from mpl_toolkits.basemap import Basemap


def test_shiftdata_rounding():
    m = Basemap(projection='cyl', llcrnrlon=-180, urcrnrlon=180,
                llcrnrlat=-90, urcrnrlat=90, lon_0=0)

    lonsin = np.array([104.123456789, -75.987654321])
    datain = np.array([1.0, 2.0])

    lonsout, dataout = m.shiftdata(lonsin, datain)

    expected_lons = np.round(lonsin, 6)

    # Compare sorted arrays to avoid order sensitivity
    np.testing.assert_array_equal(np.sort(lonsout), np.sort(expected_lons))
    assert np.array_equal(dataout, datain)
