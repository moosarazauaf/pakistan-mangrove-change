"""Where do CGMD and GMW disagree on 1996 to 2020 change, and how well was each period observed?"""
import os, sys, time
sys.path.insert(0, __file__.rsplit("scripts", 1)[0] + "src")
from common import *
import pandas as pd
ee = init()
OUT = os.path.join(__file__.rsplit("scripts", 1)[0], "outputs")
keep = pakistan_coast(ee)
area = ee.Image.pixelArea().divide(1e4)

def call(f, n=4):
    for i in range(n):
        try: return f()
        except Exception as e:
            if i == n - 1: raise
            print("retry:", str(e)[:120]); time.sleep(8 * (i + 1))

def gmw(y): return ee.ImageCollection(GMW).filterDate(f"{y}-01-01", f"{y+1}-01-01").mosaic().gt(0).unmask(0)
def cg(y): return fcc_year(ee, y).gt(0).unmask(0)
for g_end in (2020, 2019):
    code = cg(1996).multiply(8).add(cg(2020).multiply(4)).add(gmw(1996).multiply(2)).add(gmw(g_end)).toInt().rename("code").updateMask(keep)
    rows = []
    for name, prov, w, e in SITES:
        r = call(lambda: area.addBands(code).reduceRegion(ee.Reducer.sum().group(1, "code"), site_geom(ee, w, e), 30, maxPixels=1e10, tileScale=8).getInfo())
        rows += [dict(site=name, code=x["code"], ha=x["sum"]) for x in r["groups"] if x["code"] > 0]
    d = pd.DataFrame(rows); d.to_csv(f"{OUT}/transition_codes_gmw{g_end}.csv", index=False)
    t = d.groupby("code").ha.sum()
    print(f"GMW end {g_end}: total union ha {t.sum():.0f}")
    print("  CGMD 96", round(t[[c for c in t.index if c & 8]].sum()), "20", round(t[[c for c in t.index if c & 4]].sum()), "| GMW 96", round(t[[c for c in t.index if c & 2]].sum()), g_end, round(t[[c for c in t.index if c & 1]].sum()))
    for c, v in t.sort_values(ascending=False).head(10).items():
        print(f"  C96={c>>3&1} C20={c>>2&1} G96={c>>1&1} G{g_end}={c&1}: {v:.0f}")

delta = site_geom(ee, 67.2, 68.6).intersection(ee.Geometry.Rectangle([67.2, 23.5, 68.6, 25.0]), 100)
rows = []
for col, nm in [("LANDSAT/LT05/C02/T1_L2", "L5"), ("LANDSAT/LE07/C02/T1_L2", "L7"), ("LANDSAT/LC08/C02/T1_L2", "L8"), ("LANDSAT/LC09/C02/T1_L2", "L9")]:
    h = call(lambda: ee.ImageCollection(col).filterBounds(delta).map(lambda i: i.set("yr", ee.Date(i.get("system:time_start")).get("year"))).aggregate_histogram("yr").getInfo())
    rows += [dict(sensor=nm, year=int(float(k)), scenes=v) for k, v in h.items()]
s = pd.DataFrame(rows); s.to_csv(f"{OUT}/landsat_scene_counts.csv", index=False)
print(s.pivot_table(index="year", values="scenes", aggfunc="sum").T.to_string())
