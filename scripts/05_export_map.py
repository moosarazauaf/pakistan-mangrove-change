"""Export the product-disagreement classes over the Indus Delta for the published map."""
import os, sys, urllib.request
sys.path.insert(0, __file__.rsplit("scripts", 1)[0] + "src")
from common import *
ee = init()
OUT = os.path.join(__file__.rsplit("scripts", 1)[0], "outputs")
gmw = lambda y: ee.ImageCollection(GMW).filterDate(f"{y}-01-01", f"{y+1}-01-01").mosaic().gt(0).unmask(0)
cg = lambda y: fcc_year(ee, y).gt(0).unmask(0)
code = cg(1996).multiply(8).add(cg(2020).multiply(4)).add(gmw(1996).multiply(2)).add(gmw(2020)).toInt()
# 1 both maps throughout, 2 gain in CGMD, 3 GMW only throughout, 4 loss in GMW only, 5 any other mapped pixel
cls = (code.gt(0).multiply(5)
       .where(code.bitwiseAnd(12).eq(4), 2)
       .where(code.eq(3), 3)
       .where(code.eq(2), 4)
       .where(code.eq(15), 1)
       .toByte().updateMask(pakistan_coast(ee)).unmask(0).rename("class"))
url = cls.getDownloadURL({"region": ee.Geometry.Rectangle([66.95, 23.55, 68.45, 25.05], None, False), "crs": "EPSG:4326", "dimensions": "3000x3000", "format": "GEO_TIFF"})
path = f"{OUT}/disagreement_indus_delta.tif"
urllib.request.urlretrieve(url, path)
print(path, os.path.getsize(path) // 1024, "KB")
# The real shoreline and creeks for the map background: water in at least half of JRC observations.
# JRC leaves out the open sea, so anything outside the LSIB land polygons counts as water too.
land = ee.Image.constant(0).paint(ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017"), 1)
water = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").unmask(0).gte(50).Or(land.eq(0)).toByte().rename("water")
url = water.getDownloadURL({"region": ee.Geometry.Rectangle([66.95, 23.55, 68.45, 25.05], None, False), "crs": "EPSG:4326", "dimensions": "3000x3000", "format": "GEO_TIFF"})
urllib.request.urlretrieve(url, f"{OUT}/water_indus_delta.tif")
