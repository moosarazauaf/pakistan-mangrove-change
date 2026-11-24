"""Annual canopy-equivalent and mapped area per site from CGMD, and GMW v3 extent."""
import json, os, sys, time
sys.path.insert(0, __file__.rsplit("scripts", 1)[0] + "src")
from common import *
import pandas as pd
ee = init()
keep = pakistan_coast(ee)
area = ee.Image.pixelArea().divide(1e4)
OUT = os.path.join(__file__.rsplit("scripts", 1)[0], "outputs")

def call(f, n=4):
    for i in range(n):
        try: return f()
        except Exception as e:
            if i == n - 1: raise
            print("retry:", str(e)[:120]); time.sleep(8 * (i + 1))

canopy = ee.Image.cat([fcc_year(ee, y).divide(100).multiply(area).updateMask(keep).rename(f"c{y}") for y in YEARS])
mapped = ee.Image.cat([fcc_year(ee, y).gt(0).multiply(area).updateMask(fcc_year(ee, y).mask()).updateMask(keep).rename(f"m{y}") for y in YEARS])
gmw = ee.ImageCollection(GMW)
gyears = sorted(set(gmw.aggregate_array("system:time_start").map(lambda d: ee.Date(d).get("year")).getInfo()))
gimg = ee.Image.cat([gmw.filterDate(f"{y}-01-01", f"{y+1}-01-01").mosaic().gt(0).unmask(0).multiply(area).updateMask(keep).rename(f"g{y}") for y in gyears])
rows, grows = [], []
for name, prov, w, e in SITES:
    g = site_geom(ee, w, e)
    a = call(lambda: canopy.addBands(mapped).reduceRegion(ee.Reducer.sum(), g, 30, maxPixels=1e10, tileScale=8).getInfo())
    b = call(lambda: gimg.reduceRegion(ee.Reducer.sum(), g, 30, maxPixels=1e10, tileScale=8).getInfo())
    rows += [dict(site=name, province=prov, year=y, canopy_ha=a[f"c{y}"], mapped_ha=a[f"m{y}"]) for y in YEARS]
    grows += [dict(site=name, province=prov, year=y, gmw_ha=b[f"g{y}"]) for y in gyears]
    print(name, "canopy 1984", round(a["c1984"]), "2023", round(a["c2023"]), "| mapped", round(a["m1984"]), round(a["m2023"]), "| GMW", gyears[0], round(b[f"g{gyears[0]}"]), gyears[-1], round(b[f"g{gyears[-1]}"]), flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/annual_area_by_site.csv", index=False)
pd.DataFrame(grows).to_csv(f"{OUT}/gmw_area_by_site.csv", index=False)
