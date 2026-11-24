"""Shared Earth Engine definitions for the Pakistan mangrove study."""
import os
import sys

FLOOD_MAP = r"C:/Users/Moosa/OneDrive/Desktop/Flood map"
CGMD = "projects/mangrovedatahub2/assets/CGMD-AFCC30"
GMW = "projects/sat-io/open-datasets/GMW/extent/GMW_V3"
YEARS = list(range(1984, 2024))
COAST = [61.4, 23.5, 68.6, 26.0]


def init():
    sys.path.insert(0, FLOOD_MAP + "/backend")
    from dotenv import load_dotenv
    load_dotenv(FLOOD_MAP + "/backend/.env")
    load_dotenv(FLOOD_MAP + "/.env")
    cwd = os.getcwd()
    os.chdir(FLOOD_MAP + "/backend")  # the key path in .env is relative to the backend
    from gee_layers import init_ee
    init_ee()
    os.chdir(cwd)
    import ee
    return ee


# Mangrove systems of Pakistan's coast, west to east. Longitude limits sit in
# the gaps between systems in GMW 2020.
SITES = [
    ("Jiwani", "Balochistan", 61.40, 62.00),
    ("Kalmat Khor", "Balochistan", 63.70, 64.50),
    ("Miani Hor", "Balochistan", 66.20, 66.72),
    ("Karachi coast", "Sindh", 66.72, 67.20),
    ("Indus Delta", "Sindh", 67.20, 68.60),
]


def pakistan_coast(ee):
    """The coastal box with Iran and India removed (LSIB land borders)."""
    box = ee.Geometry.Rectangle(COAST)
    lsib = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017").filter(ee.Filter.inList("country_na", ["Iran", "India"]))
    return ee.Image.constant(1).clip(box).paint(lsib, 0).selfMask()


def site_geom(ee, west, east):
    return ee.Geometry.Rectangle([west, COAST[1], east, COAST[3]], None, False)


def fcc_year(ee, y):
    """FCC for one year across both CGMD tiles; the MiddleEast tile ends at 67E."""
    c = ee.ImageCollection(CGMD).filter(ee.Filter.eq("year", y)).filter(ee.Filter.inList("tile_id", ["MiddleEast", "SouthAsia"]))
    return c.select("FCC").mosaic().setDefaultProjection(c.first().select("FCC").projection())
