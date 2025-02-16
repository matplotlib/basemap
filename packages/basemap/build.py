def get_requirements(filename):
    """Read requirements from file."""
    with open(filename) as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

def get_package_data():
    return {
        "mpl_toolkits.basemap": ["data/*"]
    }