# SI disclosure: the code, the data analysis and the text of this project were produced with Claude (Anthropic), an AI assistant, under the author's direction. No SI-generated images, video or audio are used; every chart is drawn by code from NASA data.
"""Rasterise admin boundaries (geoBoundaries, CC BY 4.0) onto the 0.1 degree grid and write simplified outlines."""
import json, sys, numpy as np
from shapely.geometry import shape, box, mapping
from shapely import contains_xy
WEST, EAST, SOUTH, NORTH, RES = 94.5, 119.5, -6.0, 7.5, 0.1
NX, NY = 250, 135
win = box(WEST, SOUTH, EAST, NORTH)
xs = WEST + (np.arange(NX) + 0.5) * RES
ys = SOUTH + (np.arange(NY) + 0.5) * RES
XX, YY = np.meshgrid(xs, ys)
CTRY = {"IDN": "Indonesia", "MYS": "Malaysia", "BRN": "Brunei", "SGP": "Singapore"}
regions = []          # id -> dict
grids = {1: np.full((NY, NX), -1, np.int16), 2: np.full((NY, NX), -1, np.int16)}
feats_out = []
def fix_name(s):
    s = s.strip()
    return s.title() if s.isupper() else s
for iso in ["IDN", "MYS", "BRN", "SGP"]:
    for lvl in (1, 2):
        if iso in ("BRN", "SGP") and lvl == 2:
            continue
        g = json.load(open("geo/%s_ADM%d.geojson" % (iso, lvl)))
        for f in g["features"]:
            geom = shape(f["geometry"])
            if not geom.is_valid:
                geom = geom.buffer(0)
            if not geom.intersects(win):
                continue
            x0, y0, x1, y1 = geom.bounds
            ix0 = max(0, int((x0 - WEST) / RES)); ix1 = min(NX, int((x1 - WEST) / RES) + 1)
            iy0 = max(0, int((y0 - SOUTH) / RES)); iy1 = min(NY, int((y1 - SOUTH) / RES) + 1)
            if ix1 <= ix0 or iy1 <= iy0:
                continue
            m = contains_xy(geom, XX[iy0:iy1, ix0:ix1], YY[iy0:iy1, ix0:ix1])
            n = int(m.sum())
            if n == 0:
                continue
            # For Brunei and Singapore treat the whole country as one level-1 region
            rid = len(regions)
            sub = grids[lvl][iy0:iy1, ix0:ix1]
            sub[m & (sub < 0)] = rid
            lat_c = float(YY[iy0:iy1, ix0:ix1][m].mean())
            area = float((np.cos(np.radians(YY[iy0:iy1, ix0:ix1][m])) * (111.32 * RES) * (110.57 * RES)).sum())
            simp = geom.intersection(win).simplify(0.02 if lvl == 1 else 0.03, preserve_topology=True)
            regions.append(dict(id=rid, name=fix_name(f["properties"]["shapeName"]), country=CTRY[iso], iso=iso,
                                level=lvl, ncell=n, area_km2=round(area), parent=None,
                                lon=round(float(XX[iy0:iy1, ix0:ix1][m].mean()), 2), lat=round(lat_c, 2)))
            feats_out.append(simp)

from shapely.ops import unary_union
def dilate(g):
    """Coastal cells whose centre falls in the sea inherit the neighbouring land region (one cell only)."""
    out = g.copy()
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx == 0 and dy == 0: continue
            sh = np.full_like(g, -1)
            ys0, ys1 = max(0, dy), NY + min(0, dy); xs0, xs1 = max(0, dx), NX + min(0, dx)
            sh[ys0:ys1, xs0:xs1] = g[ys0 - dy:ys1 - dy, xs0 - dx:xs1 - dx]
            m = (out < 0) & (sh >= 0)
            out[m] = sh[m]
    return out
grids[1], grids[2] = dilate(grids[1]), dilate(grids[2])
g1, g2 = grids[1], grids[2]
for r in regions:
    if r["level"] == 2:
        p = g1[g2 == r["id"]]; p = p[p >= 0]
        r["parent"] = int(np.bincount(p).argmax()) if len(p) else None
DROP = {"Banten", "Central Java", "West Java", "East Java", "Central Sulawesi", "South Sulawesi", "West Sulawesi"}
MACRO = {
 "Sumatra": ["Aceh","Bengkulu","Jambi","Lampung","Riau","West Sumatra","South Sumatra","North Sumatra","Bangka-Belitung Islands","Riau Islands"],
 "Kalimantan": ["Central Kalimantan","South Kalimantan","West Kalimantan","East Kalimantan","North Kalimantan"],
 "Peninsular Malaysia": ["Selangor","Johor","Kuala Lumpur","Malacca","Negeri Sembilan","Kelantan","Terengganu","Pahang","Kedah","Perlis","Perak","Penang","Putrajaya"],
 "Malaysian Borneo": ["Sabah","Sarawak","Labuan"],
}
final, cells, geoms = [], [], []
def add(name, country, level, cellidx, geom, parent=None):
    lat = SOUTH + (cellidx // NX + 0.5) * RES; lon = WEST + (cellidx % NX + 0.5) * RES
    area = float((np.cos(np.radians(lat)) * (111.32 * RES) * (110.57 * RES)).sum())
    if geom is not None and not geom.is_empty:   # true polygon area (inside the window), not the dilated cell count
        area = geom.area * np.cos(np.radians(float(lat.mean()))) * 111.32 * 110.57
    final.append(dict(id=len(final), name=name, country=country, level=level, ncell=int(len(cellidx)),
                      area_km2=round(area), parent=parent, lon=round(float(lon.mean()), 2), lat=round(float(lat.mean()), 2)))
    cells.append(cellidx.astype(np.int32)); geoms.append(geom)
    return final[-1]["id"]
flat1, flat2 = g1.ravel(), g2.ravel()
old2new = {}
keep1 = [r for r in regions if r["level"] == 1 and r["name"] not in DROP and r["iso"] in ("IDN", "MYS")]
# level 0 macro regions
allc = []
for mname, members in MACRO.items():
    ids = [r["id"] for r in keep1 if r["name"] in members]
    ci = np.where(np.isin(flat1, ids))[0]; allc.append(ci)
    add(mname, "Malaysia" if "Malays" in mname else "Indonesia", 0, ci, unary_union([feats_out[i] for i in ids]).simplify(0.03))
bru = [r["id"] for r in regions if r["level"] == 1 and r["iso"] == "BRN"]
ci_b = np.where(np.isin(flat1, bru))[0]
whole = np.unique(np.concatenate(allc + [ci_b]))
bgeom = unary_union([feats_out[i] for i in bru])
add("Whole region", "Indonesia, Malaysia, Brunei", 0, whole, unary_union([g for g in geoms] + [bgeom]))
add("Borneo (whole island)", "Indonesia, Malaysia, Brunei", 0, np.unique(np.concatenate([allc[1], allc[3], ci_b])), unary_union([geoms[1], geoms[3], bgeom]))
macro_of = {m: i for i, (mname, members) in enumerate(MACRO.items()) for m in members}
for r in keep1:
    ci = np.where(flat1 == r["id"])[0]
    if len(ci) < 4: continue
    old2new[r["id"]] = add(r["name"], r["country"], 1, ci, feats_out[r["id"]], parent=macro_of.get(r["name"]))
old2new["BRN"] = add("Brunei", "Brunei", 1, ci_b, unary_union([feats_out[i] for i in bru]).simplify(0.02), parent=5)
for r in regions:
    if r["level"] == 2 and r["parent"] in old2new:
        ci = np.where(flat2 == r["id"])[0]
        if len(ci) < 4: continue
        add(r["name"], r["country"], 2, ci, feats_out[r["id"]], parent=old2new[r["parent"]])
def rnd(o):
    if isinstance(o, float): return round(o, 3)
    if isinstance(o, (list, tuple)): return [rnd(v) for v in o]
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    return o
feats = [dict(type="Feature", properties=dict(id=r["id"], name=r["name"], level=r["level"]), geometry=rnd(mapping(g)))
         for r, g in zip(final, geoms) if r["level"] >= 1]
json.dump(dict(type="FeatureCollection", features=feats), open("web/data/regions.geojson", "w"), separators=(",", ":"))
landmask = np.zeros(NX * NY, np.uint8); landmask[whole] = 1
np.save("pipeline/landmask.npy", landmask)
np.savez("pipeline/region_cells.npz", **{str(i): c for i, c in enumerate(cells)})
json.dump(final, open("pipeline/regions.json", "w"))
import os
print(len(final), "regions:", [sum(r["level"] == l for r in final) for l in (0, 1, 2)], "| land cells", int(landmask.sum()),
      "| regions.geojson MB", round(os.path.getsize("web/data/regions.geojson") / 1e6, 2))
print([ (r["name"], r["ncell"]) for r in final if r["level"] == 0])
