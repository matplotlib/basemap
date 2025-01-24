import os
import sys
from pathlib import Path


def build_geos():
    GEOS_DIR = os.environ.get("GEOS_DIR", str(Path.home() / ".local" / "geos"))

    try:
        import utils

        geos = utils.GeosLibrary("3.6.5")
        geos.build(installdir=GEOS_DIR)

        # Add GEOS directory to environment variables
        os.environ["GEOS_DIR"] = GEOS_DIR
        os.environ["LD_LIBRARY_PATH"] = (
            f"{GEOS_DIR}/lib:{os.environ.get('LD_LIBRARY_PATH', '')}"
        )

    except Exception as e:
        print(f"Error building GEOS: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    build_geos()
