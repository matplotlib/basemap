import numpy as np
import pytest
from mpl_toolkits.basemap import Basemap

# Extracted testable version of the core plotting logic
def compute_arc_points(m, lon1, lat1, lon2, lat2, arc_height_ratio=0.15):
    t = np.linspace(0, 1, 300)
    lons = lon1 + (lon2 - lon1) * t
    lats = lat1 + (lat2 - lat1) * t
    mid_lat = (lat1 + lat2) / 2
    arc_direction = -1 if mid_lat > 0 else 1
    arc_height = arc_height_ratio * abs(lat2 - lat1 + 1e-6)
    bulge = np.cos(np.pi * (t - 0.5))
    lats += arc_direction * arc_height * bulge
    lats = np.clip(lats, -85, 85)
    x, y = m(lons, lats)
    return x, y

# === Unit test ===
def test_better_great_circle():
    m = Basemap(projection='merc',
                llcrnrlon=-180, urcrnrlon=180,
                llcrnrlat=-60, urcrnrlat=80,
                lat_ts=20, resolution='c')
    # Example input: Anchorage → Moscow
    x, y = compute_arc_points(m, -149.9003, 61.2181, 37.6173, 55.7558)

    # Test: All values are finite
    assert np.all(np.isfinite(x)), "X coordinates contain non-finite values"
    assert np.all(np.isfinite(y)), "Y coordinates contain non-finite values"

    # Test: Array has expected length
    assert len(x) == 300 and len(y) == 300, "Coordinate arrays are incorrect length"

    # Test: Latitude transformation is within Mercator-safe bounds
    _, test_lats = m(x, y, inverse=True)
    assert np.all((test_lats >= -85) & (test_lats <= 85)), "Some latitudes are out of bounds"