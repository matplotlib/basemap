import itertools


def get_data_files():
    # Define some helper lists
    basenames = [
        "countries",
        "countriesmeta",
        "gshhs",
        "gshhsmeta",
        "rivers",
        "riversmeta",
        "states",
        "statesmeta",
    ]
    resolutions = ["c", "l", "i", "h", "f"]
    grids = ["1.25", "2.5", "5", "10"]

    # Define data assets
    data_dat_files = [
        f"{basename}_{res}.dat"
        for basename, res in itertools.product(basenames, resolutions[:3])
    ]
    data_bin_files = [
        f"lsmask_{grid}min_{res}.bin"
        for grid, res in itertools.product(grids, resolutions)
    ]
    data_usc_files = [f"UScounties.{ext}" for ext in ("dbf", "prj", "shp", "shx")]
    data_other_files = ["epsg", "bmng.jpg", "etopo1.jpg", "shadedrelief.jpg"]

    return data_dat_files + data_bin_files + data_usc_files + data_other_files
