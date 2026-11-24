"""The study's two figures: where the products disagree, and what each says over time.

Run: uv run --no-project --with matplotlib --with rasterio --with numpy --with pandas --with shapely python scripts/06_figures.py
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon as MPoly, Rectangle
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import shape

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)
PROV = r"C:/Users/Moosa/OneDrive/Desktop/moosarazauaf.github.io/public/geo/pakistan-provinces.geojson"
R = json.load(open(f"{OUT}/results.json"))
plt.rcParams.update({"font.family": "Segoe UI", "font.size": 11, "axes.edgecolor": "#232323", "axes.linewidth": 0.7})

INK, MUTED, GRID = "#111111", "#5c625e", "#e4e4df"
BOTH, GAIN, GMW_ONLY, LOSS = "#007350", "#35c4a0", "#7a6bd6", "#e0701a"
LAND, SEA = "#ece9e0", "#f7f9fa"
provinces = [shape(f["geometry"]) for f in json.load(open(PROV, encoding="utf8"))["features"]]


def draw_land(ax, lw=0.5):
    for g in provinces:
        for poly in g.geoms:
            ax.add_patch(MPoly(np.array(poly.exterior.coords), closed=True, fc=LAND, ec="#8f928c", lw=lw, zorder=1))


# ---- Figure 1: the map ----
with rasterio.open(f"{OUT}/disagreement_indus_delta.tif") as r:
    a, b = r.read(1), r.bounds
fig = plt.figure(figsize=(12.8, 8), dpi=200, facecolor="white")
ax = fig.add_axes([0.03, 0.06, 0.60, 0.87], facecolor=LAND, anchor="W")
with rasterio.open(f"{OUT}/water_indus_delta.tif") as r:
    w = r.read(1)
ax.imshow(np.ma.masked_equal(w, 0), extent=(b.left, b.right, b.bottom, b.top), cmap=ListedColormap(["#d6e4ec"]), interpolation="nearest", zorder=1)
show = np.ma.masked_equal(np.where(a == 5, 0, a), 0)
ax.imshow(show, extent=(b.left, b.right, b.bottom, b.top), cmap=ListedColormap([BOTH, GAIN, GMW_ONLY, LOSS]), vmin=0.5, vmax=4.5, interpolation="nearest", zorder=2)
x0, x1, y0, y1 = 66.95, 68.45, 23.62, 25.02
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect(1 / np.cos(np.radians(24.3)))
ax.xaxis.tick_top(); ax.yaxis.tick_right()
ax.set_xticks(np.arange(67.0, 68.5, 0.25)); ax.set_yticks(np.arange(23.75, 25.1, 0.25))
ax.set_xticklabels([f"{v:.2f}°E" for v in ax.get_xticks()]); ax.set_yticklabels([f"{v:.2f}°N" for v in ax.get_yticks()], rotation=-90, va="center")
ax.tick_params(direction="in", length=5, width=0.7, labelsize=10, colors=INK)
for name, lon, lat in [("Karachi", 67.03, 24.9), ("Keti Bandar", 67.45, 24.15), ("Shah Bandar", 67.92, 24.17)]:
    ax.plot(lon, lat, "o", ms=4, mfc="white", mec=INK, mew=0.8, zorder=4)
    ax.annotate(name, (lon, lat), xytext=(6, 5), textcoords="offset points", fontsize=10.5, color=INK, zorder=5,
                path_effects=[matplotlib.patheffects.withStroke(linewidth=2.6, foreground="white")])
ax.text(67.0, 23.86, "Arabian Sea", fontsize=13, style="italic", color="#54707f", zorder=3)

d = R["disputed"]
hk = lambda k: f"{d[k]['ha'] / 1000:.1f}k ha"
gain_ha = sum(v for c, v in pd.read_csv(f"{OUT}/transition_codes_gmw2020.csv").groupby("code").ha.sum().items() if (c & 12) == 4) / 1000
leg = fig.legend(handles=[
    Patch(color=BOTH, label=f"In both maps, 1996 and 2020 ({hk('both_keep')})"),
    Patch(color=GAIN, label=f"New in CGMD since 1996 ({gain_ha:.1f}k ha)"),
    Patch(color=GMW_ONLY, label=f"In GMW only, both dates ({hk('gmw_only_stable')})"),
    Patch(color=LOSS, label=f"Lost in GMW, never in CGMD ({hk('gmw_loss_cgmd_absent')})"),
], loc="upper left", bbox_to_anchor=(0.658, 0.43), frameon=False, fontsize=11.5, title="What two global maps say, 1996 to 2020", title_fontsize=12, alignment="left", labelspacing=0.8)

km = 25
deg = km / (111.32 * np.cos(np.radians(24.3)))
sx, sy = 67.0, 23.70
ax.plot([sx, sx + deg, sx + 2 * deg], [sy] * 3, color=INK, lw=1.1, zorder=6)
for i in range(3):
    ax.plot([sx + i * deg] * 2, [sy, sy + 0.012], color=INK, lw=1.1, zorder=6)
    ax.text(sx + i * deg, sy + 0.018, f"{i * km}" + (" km" if i == 2 else ""), ha="center", fontsize=9.5, zorder=6)
ax.annotate("N", xy=(68.36, 24.96), xytext=(68.36, 24.84), ha="center", fontsize=12, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2), zorder=6)

ins = fig.add_axes([0.665, 0.49, 0.30, 0.44], facecolor=SEA)
draw_land(ins, lw=0.35)
ins.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec=LOSS, lw=1.3, zorder=4))
for lon, lat in [(61.75, 25.15), (64.1, 25.4), (66.45, 25.45), (66.95, 24.85), (67.6, 24.2)]:
    ins.plot(lon, lat, "o", ms=3.2, color=BOTH, mec="white", mew=0.4, zorder=5)
ins.set_xlim(60.5, 78); ins.set_ylim(23, 37.3); ins.set_aspect(1 / np.cos(np.radians(30))); ins.set_xticks([]); ins.set_yticks([])
fig.text(0.665, 0.455, "Green dots: the five mangrove systems studied. Box: the main map.", fontsize=9.5, color=MUTED)
fig.text(0.03, 0.018, "CGMD (Zhang et al. 2026), Global Mangrove Watch v3 (Bunting et al. 2022), water from JRC Global Surface Water. WGS 84. Boundaries: geoBoundaries; no position on disputed territory is implied.", fontsize=8.5, color=MUTED)
fig.savefig(f"{FIG}/fig1_disagreement_map.png", facecolor="white")
plt.close(fig)

# ---- Figure 2: the series ----
s = pd.read_csv(f"{OUT}/annual_area_by_site.csv").groupby("year").mapped_ha.sum() / 1000
g = pd.Series({int(k): v / 1000 for k, v in R["gmw"].items()})
ex = json.load(open(f"{OUT}/extent_checks.json"))
A = R["arbiter"]
fig, ax = plt.subplots(figsize=(12.8, 8), dpi=200, facecolor="white")
fig.subplots_adjust(left=0.07, right=0.83, top=0.95, bottom=0.08)
ax.yaxis.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.plot(s.index, s.values, color=BOTH, lw=2, zorder=3)
ax.plot(g.index, g.values, color=GMW_ONLY, lw=2, ls=(0, (4, 2)), marker="o", ms=6, mec="white", mew=1, zorder=3)
yrs = [1996, 2020]
vals = np.array([A["ha_1996_epoch"], A["ha_2020_epoch"]]) / 1000
errs = np.array([[A["ha_1996_epoch"] - A["ci_1996"][0], A["ha_2020_epoch"] - A["ci_2020"][0]], [A["ci_1996"][1] - A["ha_1996_epoch"], A["ci_2020"][1] - A["ha_2020_epoch"]]]) / 1000
ax.plot(yrs, vals, color=INK, lw=1, ls=":", zorder=4)
ax.errorbar(yrs, vals, yerr=errs, fmt="s", ms=9, color=INK, mec="white", mew=1.2, capsize=4, lw=1.3, zorder=5)
ax.plot([2000], [ex["giri_2000"]["ha"] / 1000], "D", ms=8, color=LOSS, mec="white", mew=1, zorder=5)
ax.plot([2020], [ex["worldcover_2020"]["ha"] / 1000], "^", ms=10, color=LOSS, mec="white", mew=1, zorder=5)
ax.annotate("Giri et al. 2000", (2000, ex["giri_2000"]["ha"] / 1000), xytext=(8, -14), textcoords="offset points", fontsize=10.5, color=INK)
ax.annotate("ESA WorldCover 2020", (2020, ex["worldcover_2020"]["ha"] / 1000), xytext=(-8, 9), textcoords="offset points", ha="right", fontsize=10.5, color=INK)
ax.text(1996.9, 63.4, f"Landsat check: {A['change_pct']:+.0f}%, 1996 to 2020", fontsize=11, color=INK)
gp = 100 * (g[2020] - g[1996]) / g[1996]
cp = 100 * (R["national"]["mapped_2020"] - R["national"]["mapped_1996"]) / R["national"]["mapped_1996"]
ax.text(2023.5, s[2023], f"CGMD\n{cp:+.0f}%, 1996 to 2020", color=INK, va="center", fontsize=11)
ax.text(2020.7, g[2020] + 3.4, f"GMW\n{gp:+.0f}%", color=INK, va="center", fontsize=11, clip_on=False)
ax.set_xlim(1983, 2023.3); ax.set_ylim(30, 110)
ax.legend(handles=[Line2D([], [], color=BOTH, lw=2, label="CGMD, annual"), Line2D([], [], color=GMW_ONLY, lw=2, ls=(0, (4, 2)), marker="o", ms=6, label="Global Mangrove Watch"),
                   Line2D([], [], color=INK, ls="", marker="s", ms=8, label="Landsat check, this study (95% interval)"), Line2D([], [], color=LOSS, ls="", marker="D", ms=7, label="Other single-year maps")],
          loc="upper left", frameon=False, fontsize=10.5)
ax.set_ylabel("Mangrove area in Pakistan (thousand ha)", color=INK)
ax.tick_params(length=3, width=0.7, colors=INK)
fig.savefig(f"{FIG}/fig2_area_series.png", facecolor="white")
plt.close(fig)
print("figures written")
