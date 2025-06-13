import pytest
from datetime import datetime
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


def test_nightshade_runs_without_error():
    # Setup Basemap with correct lon_0 to avoid wrapping issues
    m = Basemap(projection='robin', lon_0=0, resolution='c')

    # Draw basemap
    m.drawcoastlines()
    m.fillcontinents(color='coral', lake_color='aqua')
    m.drawmapboundary(fill_color='aqua')

    # Define test datetime
    test_date = datetime(2024, 10, 10, 6)

    try:
        # Attempt to draw nightshade
        shade = m.nightshade(test_date, alpha=0.2)
    except Exception as e:
        pytest.fail(f"nightshade() raised an error: {e}")

    # Optional: Check object type
    assert hasattr(shade, "collections"), "nightshade() did not return a valid contour object"

    plt.close()  # Clean up the figure after test
