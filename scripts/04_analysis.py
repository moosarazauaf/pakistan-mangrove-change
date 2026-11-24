"""Which global product is right about Pakistan's mangroves? All reported statistics are written to outputs/results.json.

Run: uv run --no-project --with pandas --with numpy --with scipy --with pymannkendall python scripts/04_analysis.py
"""
import json, os
import numpy as np
import pandas as pd
import pymannkendall as mk

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")
rng = np.random.default_rng(42)
R = {}

# IPCC 2013 Wetlands Supplement, Chapter 4: Tables 4.2, 4.3 and 4.5, tropical dry mangrove.
AGB, AGB_LO, AGB_HI = 92, 88, 97
ROOT, ROOT_LO, ROOT_HI = 0.29, 0.28, 0.30
CF, CF_LO, CF_HI = 0.451, 0.429, 0.471
C_HA = AGB * (1 + ROOT) * CF
C_LO, C_HI = AGB_LO * (1 + ROOT_LO) * CF_LO, AGB_HI * (1 + ROOT_HI) * CF_HI
R["biomass_c_t_per_ha"] = {"central": round(C_HA, 1), "low": round(C_LO, 1), "high": round(C_HI, 1)}

# ---- 1. Trends, the tutorial's way and with autocorrelation handled ----
a = pd.read_csv(f"{OUT}/annual_area_by_site.csv")
nat = a.groupby("year")[["canopy_ha", "mapped_ha"]].sum()
nat["mean_fcc"] = 100 * nat.canopy_ha / nat.mapped_ha


def trend(s):
    o, h = mk.original_test(s.values), mk.hamed_rao_modification_test(s.values)
    return {"sen_per_yr": round(float(o.slope), 3), "p_original": float(f"{o.p:.2e}"), "p_hamed_rao": float(f"{h.p:.2e}"), "tau": round(float(o.Tau), 3)}


R["trends"] = {}
for label, frame in [("1984-2023", nat), ("2000-2023", nat.loc[2000:])]:
    R["trends"][label] = {c: trend(frame[c]) for c in ("canopy_ha", "mapped_ha", "mean_fcc")}
R["sites_canopy"] = {}
for site, g in a.groupby("site"):
    g = g.set_index("year")
    R["sites_canopy"][site] = {"province": g.province.iloc[0], "ha_1996": round(g.mapped_ha[1996]), "ha_2020": round(g.mapped_ha[2020]), "canopy_1984": round(g.canopy_ha[1984]), "canopy_2023": round(g.canopy_ha[2023]), **trend(g.loc[2000:].canopy_ha)}
R["national"] = {"canopy_1984": round(nat.canopy_ha[1984]), "canopy_2023": round(nat.canopy_ha[2023]), "mapped_1996": round(nat.mapped_ha[1996]), "mapped_2020": round(nat.mapped_ha[2020]), "mean_fcc_1984": round(nat.mean_fcc[1984], 1), "mean_fcc_2023": round(nat.mean_fcc[2023], 1)}

gm = pd.read_csv(f"{OUT}/gmw_area_by_site.csv")
gn = gm.groupby("year").gmw_ha.sum()
R["gmw"] = {str(y): round(v) for y, v in gn.items()}
R["gmw_sites"] = {s: {"ha_1996": round(g.set_index("year").gmw_ha[1996]), "ha_2020": round(g.set_index("year").gmw_ha[2020])} for s, g in gm.groupby("site")}

# Does CGMD's area move with how many Landsat scenes existed that year?
sc = pd.read_csv(f"{OUT}/landsat_scene_counts.csv").groupby("year").scenes.sum()
j = pd.concat([nat.canopy_ha, sc], axis=1).dropna()
R["scene_count_check"] = {"years": len(j), "r_levels": round(float(j.corr().iloc[0, 1]), 3), "r_year_to_year_change": round(float(j.diff().dropna().corr().iloc[0, 1]), 3)}

# ---- 2. Landsat NDVI as an arbiter on the pixels the products dispute ----
p = pd.read_csv(f"{OUT}/sample_points.csv").dropna(subset=["ndvi_a", "ndvi_b"])
p["code"] = p.code.astype(int)
for bit, nm in [(3, "c96"), (2, "c20"), (1, "g96"), (0, "g20")]:
    p[nm] = np.where(p.code == 16, 0, (p.code.values >> bit) & 1)
W = pd.read_csv(f"{OUT}/transition_codes_gmw2020.csv").set_index(["site", "code"]).ha
R["n_points"] = int(len(p))
R["n_mapped_points"] = int((p.code < 16).sum())


def thresholds(pos, neg):
    grid = np.arange(-0.1, 0.8, 0.002)
    tpr = np.array([(pos >= t).mean() for t in grid])
    fpr = np.array([(neg >= t).mean() for t in grid])
    youden = grid[np.argmax(tpr - fpr)]
    # rank-sum AUC
    x = np.r_[pos, neg]
    r = pd.Series(x).rank().values[: len(pos)]
    auc = (r.sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))
    return {"youden": round(float(youden), 3), "high": round(float(np.quantile(pos, 0.05)), 3), "low": round(float(np.quantile(neg, 0.95)), 3), "auc": round(float(auc), 3), "tpr": round(float((pos >= youden).mean()), 3), "fpr": round(float((neg >= youden).mean()), 3)}


neg = p[p.code == 16]
T = {"a": thresholds(p[(p.c96 == 1) & (p.g96 == 1)].ndvi_a.values, neg.ndvi_a.values), "b": thresholds(p[(p.c20 == 1) & (p.g20 == 1)].ndvi_b.values, neg.ndvi_b.values)}
R["ndvi_thresholds"] = T
m = p[p.code < 16].copy()
keys = [k for k in W.index if k in set(zip(m.site, m.code))]
R["mapped_ha_without_sample"] = round(float(W.sum() - W.loc[keys].sum()), 1)
groups = {k: g for k, g in m.groupby(["site", "code"])}


def estimate(ta, tb, boot=0, need_wc=False, sites=None):
    ks = [k for k in groups if sites is None or k[0] in sites]
    va = [groups[k].ndvi_a.values >= ta for k in ks]
    vb = [(groups[k].ndvi_b.values >= tb) & (groups[k].wc20.values == 1 if need_wc else True) for k in ks]
    w = np.array([W[k] for k in ks])
    point = lambda A, B: (float(np.sum(w * np.array([x.mean() for x in A]))), float(np.sum(w * np.array([x.mean() for x in B]))))
    ea, eb = point(va, vb)
    if not boot:
        return ea, eb
    out = []
    for _ in range(boot):
        ia = [rng.integers(0, len(x), len(x)) for x in va]
        out.append(point([x[i] for x, i in zip(va, ia)], [x[i] for x, i in zip(vb, ia)]))
    o = np.array(out)
    ci = lambda v: [round(float(q)) for q in np.quantile(v, [0.025, 0.975])]
    return ea, eb, ci(o[:, 0]), ci(o[:, 1]), ci(o[:, 1] - o[:, 0])


ea, eb, ca, cb, cd = estimate(T["a"]["youden"], T["b"]["youden"], boot=2000)
R["arbiter"] = {"ha_1996_epoch": round(ea), "ci_1996": ca, "ha_2020_epoch": round(eb), "ci_2020": cb, "change_ha": round(eb - ea), "ci_change": cd, "change_pct": round(100 * (eb - ea) / ea, 1)}
R["arbiter_sensitivity"] = {}
for nm in ("high", "low"):
    x, y = estimate(T["a"][nm], T["b"][nm])
    R["arbiter_sensitivity"][nm + "_threshold"] = {"ha_1996_epoch": round(x), "ha_2020_epoch": round(y), "change_ha": round(y - x), "change_pct": round(100 * (y - x) / x, 1)}
# Harshest reading: a 2020 pixel counts only if WorldCover (Sentinel-1 and 2, no Landsat) also calls it mangrove. No such filter exists for 1996, so this leans against gain.
x, y, _, _, cc = estimate(T["a"]["youden"], T["b"]["youden"], boot=2000, need_wc=True)
R["arbiter_sensitivity"]["worldcover_confirmed_2020"] = {"ha_1996_epoch": round(x), "ha_2020_epoch": round(y), "change_ha": round(y - x), "ci_change": cc, "change_pct": round(100 * (y - x) / x, 1)}
R["arbiter_by_site"] = {}
for s in sorted(m.site.unique()):
    x, y = estimate(T["a"]["youden"], T["b"]["youden"], sites={s})
    R["arbiter_by_site"][s] = {"ha_1996_epoch": round(x), "ha_2020_epoch": round(y), "change_pct": round(100 * (y - x) / x, 1) if x else None}
for prov, ss in {"Sindh": {"Karachi coast", "Indus Delta"}, "Balochistan": {"Jiwani", "Kalmat Khor", "Miani Hor"}}.items():
    x, y = estimate(T["a"]["youden"], T["b"]["youden"], sites=ss)
    R["arbiter_by_site"][prov] = {"ha_1996_epoch": round(x), "ha_2020_epoch": round(y), "change_pct": round(100 * (y - x) / x, 1)}

# The disputed classes: what fraction was vegetated in each epoch?
CLASSES = {
    "both_keep": ((1, 1, 1, 1), "Both products: mangrove throughout"),
    "cgmd_gain_gmw_stable": ((0, 1, 1, 1), "CGMD gain; GMW says it was already there"),
    "cgmd_gain_gmw_absent": ((0, 1, 0, 0), "CGMD gain; GMW never maps it"),
    "gmw_loss_cgmd_absent": ((0, 0, 1, 0), "GMW loss; CGMD never maps it"),
    "gmw_only_stable": ((0, 0, 1, 1), "GMW mangrove throughout; CGMD never maps it"),
}
R["disputed"] = {}
for k, ((c96, c20, g96, g20), label) in CLASSES.items():
    code = c96 * 8 + c20 * 4 + g96 * 2 + g20
    g = m[m.code == code]
    wt = g.groupby("site").code.transform("size")
    w = g.site.map(lambda s: W.get((s, code), 0)) / wt
    wm = lambda v: float((v * w).sum() / w.sum())
    R["disputed"][k] = {"label": label, "ha": round(float(W.xs(code, level="code").sum())), "n": int(len(g)),
                        "veg_1996_epoch": round(wm(g.ndvi_a >= T["a"]["youden"]), 3), "veg_2020_epoch": round(wm(g.ndvi_b >= T["b"]["youden"]), 3),
                        "worldcover_2020": round(wm(g.wc20), 3), "giri_2000": round(wm(g.giri), 3)}


# Agreement of each 2020 map with ESA WorldCover and of each ~2000 map with Giri, weighted by stratum area.
def share(mask_col, ref):
    g = m[m[mask_col] == 1]
    w = pd.Series([W[(s, c)] for s, c in zip(g.site, g.code)], index=g.index) / g.groupby(["site", "code"]).code.transform("size")
    return round(float((g[ref] * w).sum() / w.sum()), 3)


R["worldcover_agreement_2020"] = {"share_of_cgmd_2020_also_worldcover": share("c20", "wc20"), "share_of_gmw_2020_also_worldcover": share("g20", "wc20")}
R["giri_agreement"] = {"share_of_cgmd_1996_also_giri_2000": share("c96", "giri"), "share_of_gmw_1996_also_giri_2000": share("g96", "giri")}

# ---- 3. What the choice of product does to a carbon account ----
g96, g20, c96, c20 = gn[1996], gn[2020], nat.mapped_ha[1996], nat.mapped_ha[2020]
mt = lambda ha, d=C_HA: round(ha * d / 1e6, 2)
R["carbon"] = {"gmw_change_Mt_C": mt(g20 - g96), "cgmd_change_Mt_C": mt(c20 - c96), "arbiter_change_Mt_C": mt(eb - ea),
               "arbiter_change_range_Mt_C": [mt(eb - ea, C_LO), mt(eb - ea, C_HI)], "gap_between_products_Mt_C": mt((c20 - c96) - (g20 - g96))}

# ---- 4. Do big Indus floods show up in the delta's canopy the following year? ----
piv = a.pivot(index="year", columns="site", values="canopy_ha")
delta = piv["Indus Delta"].pct_change().dropna() * 100
ctrl = piv[["Jiwani", "Kalmat Khor", "Miani Hor"]].sum(axis=1).pct_change().dropna() * 100
after = [2011, 2023]
diff = (delta - ctrl).dropna()
obs = float(diff.loc[after].mean())
perm = np.array([rng.choice(diff.values, 2, replace=False).mean() for _ in range(20000)])
R["flood_pulse"] = {"years_after_flood": after, "delta_pct": {str(y): round(float(delta[y]), 2) for y in after}, "delta_minus_balochistan_pct": round(obs, 2),
                    "other_years_mean_pct": round(float(diff.drop(after).mean()), 2), "permutation_p_two_sided": round(float((np.abs(perm - diff.mean()) >= abs(obs - diff.mean())).mean()), 3)}

json.dump(R, open(f"{OUT}/results.json", "w"), indent=1)
print(json.dumps({k: v for k, v in R.items() if k not in ("sites_canopy",)}, indent=1))
