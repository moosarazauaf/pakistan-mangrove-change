"""How much mangrove lies east of the MiddleEast tile, and what two third-party maps say."""
import json, os, sys
sys.path.insert(0, __file__.rsplit("scripts", 1)[0] + "src")
from common import *
ee = init()
OUT = os.path.join(__file__.rsplit("scripts", 1)[0], "outputs")
keep = pakistan_coast(ee)
area = ee.Image.pixelArea().divide(1e4)
box = ee.Geometry.Rectangle(COAST)
east = ee.Image.pixelLonLat().select("longitude").gte(67)


def ha(mask):
    m = mask.unmask(0).gt(0).updateMask(keep)
    img = area.updateMask(m).rename("all").addBands(area.updateMask(m.And(east)).rename("east"))
    r = img.reduceRegion(ee.Reducer.sum(), box, 30, maxPixels=1e10, tileScale=8).getInfo()
    return {"ha": round(r["all"]), "east_of_67E_ha": round(r["east"]), "east_share": round(r["east"] / r["all"], 4)}


gmw = lambda y: ee.ImageCollection(GMW).filterDate(f"{y}-01-01", f"{y+1}-01-01").mosaic()
out = {
    "cgmd_2000": ha(fcc_year(ee, 2000)),
    "cgmd_2020": ha(fcc_year(ee, 2020)),
    "cgmd_2023": ha(fcc_year(ee, 2023)),
    "gmw_1996": ha(gmw(1996)),
    "gmw_2020": ha(gmw(2020)),
    "giri_2000": ha(ee.ImageCollection("LANDSAT/MANGROVE_FORESTS").mosaic()),
    "worldcover_2020": ha(ee.ImageCollection("ESA/WorldCover/v100").first().eq(95).selfMask()),
}
json.dump(out, open(f"{OUT}/extent_checks.json", "w"), indent=1)
for k, v in out.items():
    print(k, v)
