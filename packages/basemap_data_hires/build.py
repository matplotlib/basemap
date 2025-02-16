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

    # Define data assets
    return [
        f"{basename}_{res}.dat"
        for basename, res in itertools.product(basenames, resolutions[3:])
    ]
