import json, sys
sys.path.insert(0, __file__.rsplit("scripts", 1)[0] + "src")
from common import *
ee = init()
coast = ee.Geometry.Rectangle(COAST)
c = ee.ImageCollection(CGMD)
for tile in ("MiddleEast", "SouthAsia"):
    t = c.filter(ee.Filter.eq("tile_id", tile)).filter(ee.Filter.eq("year", 2020)).first()
    g = t.geometry().bounds().coordinates().get(0).getInfo()
    xs = [p[0] for p in g]; ys = [p[1] for p in g]
    print(tile, "bounds", round(min(xs), 3), round(max(xs), 3), round(min(ys), 3), round(max(ys), 3))
    print("  hist", json.dumps(t.select("FCC").reduceRegion(ee.Reducer.fixedHistogram(0, 110, 11), coast, 300, maxPixels=1e9, tileScale=4).getInfo())[:330])
    print("  min max", t.select("FCC").reduceRegion(ee.Reducer.minMax(), coast, 300, maxPixels=1e9, tileScale=4).getInfo())
gmw20 = ee.ImageCollection(GMW).filterDate("2020-01-01", "2021-01-01").mosaic().gt(0)
print("GMW bands", ee.ImageCollection(GMW).first().bandNames().getInfo(), ee.ImageCollection(GMW).size().getInfo())
area = ee.Image.pixelArea().divide(1e4)
lon = ee.Image.pixelLonLat().select("longitude").subtract(61).divide(0.25).floor().toInt()
out = area.updateMask(gmw20.unmask(0)).addBands(lon).reduceRegion(ee.Reducer.sum().group(1, "bin"), coast, 30, maxPixels=1e10, tileScale=8).getInfo()
print("GMW 2020 ha by 0.25 deg lon:")
tot = 0
for g in out["groups"]:
    tot += g["sum"]; print(f"  {61 + g['bin'] * .25:.2f}E {g['sum']:.0f}")
print("total", round(tot))
