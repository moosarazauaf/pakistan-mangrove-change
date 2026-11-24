"""Stratified points in every product-transition class, with Landsat NDVI for both epochs and two third-party maps."""
import os, sys, time
sys.path.insert(0, __file__.rsplit("scripts", 1)[0] + "src")
from common import *
import pandas as pd
ee = init()
OUT = os.path.join(__file__.rsplit("scripts", 1)[0], "outputs")
keep = pakistan_coast(ee)
box = ee.Geometry.Rectangle(COAST)
PER_CLASS, BATCH = 150, 400

def call(f, n=5):
    for i in range(n):
        try: return f()
        except Exception as e:
            if i == n - 1: raise
            print("retry:", str(e)[:120], flush=True); time.sleep(10 * (i + 1))

def gmw(y): return ee.ImageCollection(GMW).filterDate(f"{y}-01-01", f"{y+1}-01-01").mosaic().gt(0).unmask(0)
def cg(y): return fcc_year(ee, y).gt(0).unmask(0)

code = cg(1996).multiply(8).add(cg(2020).multiply(4)).add(gmw(1996).multiply(2)).add(gmw(2020)).toInt()
union = code.gt(0)
# Non-mangrove reference: unmapped pixels within about 1 km (33 pixels) of any mapped mangrove.
near = union.fastDistanceTransform(40).sqrt().lte(33)
strat = code.where(code.eq(0).And(near), 16).selfMask().rename("code").updateMask(keep)

F = []
for s, (name, prov, w, e) in enumerate(SITES):
    r = call(lambda: strat.stratifiedSample(numPoints=PER_CLASS, classBand="code", region=site_geom(ee, w, e), scale=30, seed=42 + s, geometries=True, tileScale=8).getInfo())
    for f in r["features"]:
        f["properties"]["site"] = name
    F += r["features"]
    print(name, "points", len(r["features"]), flush=True)

def ndvi(col, red, nir):
    def f(i):
        qa = i.select("QA_PIXEL")
        clear = qa.bitwiseAnd(0b11110).eq(0).And(i.select("QA_RADSAT").eq(0))
        s = i.select([red, nir], ["r", "n"]).multiply(0.0000275).add(-0.2)
        return s.normalizedDifference(["n", "r"]).rename("ndvi").updateMask(clear.And(s.select("r").gt(0)).And(s.select("n").gt(0)))
    return col.filterBounds(box).filter(ee.Filter.calendarRange(10, 5, "month")).map(f)

early = ndvi(ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterDate("1994-01-01", "1999-01-01"), "SR_B3", "SR_B4")
late = ndvi(ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").merge(ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")).filterDate("2018-01-01", "2023-01-01"), "SR_B4", "SR_B5")
# The 90th percentile leans on low-tide views, when the canopy is least hidden by water.
stack = ee.Image.cat([
    early.reduce(ee.Reducer.percentile([90])).rename("ndvi_a"), early.count().rename("n_a"),
    late.reduce(ee.Reducer.percentile([90])).rename("ndvi_b"), late.count().rename("n_b"),
    ee.ImageCollection("LANDSAT/MANGROVE_FORESTS").mosaic().gt(0).unmask(0).rename("giri"),
    ee.ImageCollection("ESA/WorldCover/v100").first().eq(95).unmask(0).rename("wc20"),
])
rows = []
for k in range(0, len(F), BATCH):
    fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Point(f["geometry"]["coordinates"]), {"code": f["properties"]["code"], "site": f["properties"]["site"], "pid": k + j}) for j, f in enumerate(F[k:k + BATCH])])
    out = call(lambda: stack.reduceRegions(fc, ee.Reducer.first(), 30, tileScale=8).getInfo())
    for f, src in zip(sorted(out["features"], key=lambda x: x["properties"]["pid"]), F[k:k + BATCH]):
        rows.append({**f["properties"], "lon": src["geometry"]["coordinates"][0], "lat": src["geometry"]["coordinates"][1]})
    pd.DataFrame(rows).to_csv(f"{OUT}/sample_points.csv", index=False)
    print("done", len(rows), flush=True)
