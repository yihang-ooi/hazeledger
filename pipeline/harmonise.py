"""
harmonise.py - turn per-satellite daily grids into one consistent fire record.

Input : folder with agg_T/A/S/J.csv.gz + manifest.json (from firms_fetch.py)
Output: web/data/{meta.json, series_main.bin, series_l2.bin, grid_annual.bin} and pipeline/report.json

Unit of the harmonised record: "VIIRS-equivalent detections" = what Suomi NPP VIIRS (375 m) reports.
  2012-01-20 onward : Suomi NPP VIIRS itself (observed)
  before that       : MODIS (Terra+Aqua) x calibration factor learnt on the 2012+ overlap
  Terra-only months : Terra x its own factor
"""
import sys, json, os, datetime as dt
import numpy as np, pandas as pd
from scipy import sparse, stats

DATA = sys.argv[1] if len(sys.argv) > 1 else "hazeledger_data"
RECENT = sys.argv[2] if len(sys.argv) > 2 else None
HAS_RECENT = []
OUT = "web/data"
NX, NY = 250, 135
NCELL = NX * NY
EPOCH = dt.date(2000, 1, 1)
ALPHA = 20.0            # shrinkage pseudo-count (MODIS detections)
os.makedirs(OUT, exist_ok=True)


def d2s(i):
    return str(EPOCH + dt.timedelta(days=int(i)))


# ------------------------------------------------------------------ load
land = np.load("pipeline/landmask.npy").astype(bool)
regions = json.load(open("pipeline/regions.json"))
rc = np.load("pipeline/region_cells.npz")
NR = len(regions)


def load(code):
    p = os.path.join(DATA, "agg_%s.csv.gz" % code)
    if not os.path.exists(p):
        return None
    df = pd.read_csv(p)
    # recent months from the FIRMS API: use only days after the archive's final calendar year
    for rd in (RECENT, DATA):
        q = os.path.join(rd, "recent_%s.csv.gz" % code) if rd else None
        if q and os.path.exists(q):
            last_year = (EPOCH + dt.timedelta(days=int(df.day.max()))).year
            cut = (dt.date(last_year, 12, 31) - EPOCH).days
            rec = pd.read_csv(q)
            df = pd.concat([df, rec[rec.day > cut]], ignore_index=True)
            HAS_RECENT.append(code)
            break
    n_all = int(df.nD.sum() + df.nN.sum())
    df = df[land[df.cell.values]]
    df = df.groupby(["day", "cell"], as_index=False)[["nD", "nN", "frpD", "frpN"]].sum()
    df.attrs["n_all"] = n_all
    df.attrs["n_land"] = int(df.nD.sum() + df.nN.sum())
    return df


raw = {c: load(c) for c in "TASJK"}
raw = {c: v for c, v in raw.items() if v is not None and len(v)}
for c, v in raw.items():
    print("%s: %9d detections on land (%.1f%% of window)  %s .. %s" % (
        c, v.attrs["n_land"], 100 * v.attrs["n_land"] / max(1, v.attrs["n_all"]), d2s(v.day.min()), d2s(v.day.max())))

D = int(max(v.day.max() for v in raw.values())) + 1
S0 = int(raw["S"].day.min())                 # first Suomi NPP day
A0 = int(raw["A"].day.min())                 # first Aqua day
S_END = int(raw["S"].day.max())
M_END = int(min(raw["T"].day.max(), raw["A"].day.max()))
dates = pd.to_datetime([EPOCH + dt.timedelta(days=i) for i in range(D)])
year_of = dates.year.values
month_of = dates.month.values
SEASON4 = np.array([0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3])[month_of]   # JFM, AMJ, JASO, ND


def cube(dfs, cols):
    """outer-join several sensors on (day, cell)."""
    out = None
    for name, (df, c) in dfs.items():
        t = df[["day", "cell"] + c].rename(columns={k: "%s_%s" % (name, k) for k in c})
        out = t if out is None else out.merge(t, on=["day", "cell"], how="outer")
    return out.fillna(0)


# MODIS combined (Terra + Aqua) and Terra alone
M = cube({"T": (raw["T"], ["nD", "nN"]), "A": (raw["A"], ["nD", "nN"])}, None)
M["mD"] = M.T_nD + M.A_nD
M["mN"] = M.T_nN + M.A_nN
ov = M.merge(raw["S"][["day", "cell", "nD", "nN"]].rename(columns={"nD": "sD", "nN": "sN"}),
             on=["day", "cell"], how="outer").fillna(0)
ov = ov[(ov.day >= S0) & (ov.day <= min(S_END, M_END))].copy()
ov["year"] = year_of[ov.day.values]
ov["month"] = month_of[ov.day.values]
ov["m"] = ov.mD + ov.mN
ov["s"] = ov.sD + ov.sN
ov["t"] = ov.T_nD + ov.T_nN
ov["tD"] = ov.T_nD
ov["tN"] = ov.T_nN
ix, iy = ov.cell.values % NX, ov.cell.values // NX
OV_YEARS = sorted(int(y) for y in ov.year.unique())

# region lookups: level-1 province of each cell (for validation scale)
l1_of_cell = np.full(NCELL, -1, np.int32)
for r in regions:
    if r["level"] == 1:
        l1_of_cell[rc[str(r["id"])]] = r["id"]
ov["l1"] = l1_of_cell[ov.cell.values]


# ------------------------------------------------------------------ calibration models
def zone_ids(zdeg):
    if zdeg is None:
        return np.zeros(len(ov), np.int32), 1
    n = int(round(zdeg / 0.1))
    zx, zy = ix // n, iy // n
    nzx = (NX + n - 1) // n
    return (zy * nzx + zx).astype(np.int32), nzx * ((NY + n - 1) // n)


INT_EDGES = np.array([5.0, 20.0, 100.0])


def intensity_class(cell, day, x, scale=1.0):
    """How hard was the surrounding 0.5 degree block burning that month, judged from MODIS alone? 4 classes."""
    cx, cy = cell % NX, cell // NX
    z5 = (cy // 5) * 50 + cx // 5
    ym = year_of[day].astype(np.int64) * 12 + month_of[day]
    key = z5.astype(np.int64) * 100000 + ym
    u, inv = np.unique(key, return_inverse=True)
    tot = np.bincount(inv, weights=x, minlength=len(u))[inv]
    return np.searchsorted(INT_EDGES * scale, tot, side="right").astype(np.int32)


def season_of(kind, cell, day, x, scale=1.0):
    if kind == "none":
        return np.zeros(len(day), np.int32), 1
    if kind == "s4":
        return SEASON4[day].astype(np.int32), 4
    if kind == "m12":
        return (month_of[day] - 1).astype(np.int32), 12
    if kind == "int4":
        return intensity_class(cell, day, x, scale), 4
    if kind == "s4int":
        return (SEASON4[day] * 4 + intensity_class(cell, day, x, scale)).astype(np.int32), 16
    raise ValueError(kind)


class Model:
    """Ratio-of-sums calibration with shrinkage toward the regional factor.
       structure 'tot': H = k * (mD + mN);  'dn': H = kD * mD + kN * mN"""

    def __init__(self, zdeg, season, structure, pred=("mD", "mN")):
        self.zdeg, self.season, self.structure, self.pred = zdeg, season, structure, pred
        self.z, self.nz = zone_ids(zdeg)
        self.scale = 1.0 if pred[0] == "mD" else float(ov[pred[0]].sum() + ov[pred[1]].sum()) / float(ov.m.sum())
        self.s, self.ns = season_of(season, ov.cell.values, ov.day.values, ov[pred[0]].values + ov[pred[1]].values, self.scale)
        self.g = self.z * self.ns + self.s
        self.ng = self.nz * self.ns
        yi = ov.year.values - OV_YEARS[0]
        ny = len(OV_YEARS)
        self.parts = {}
        pairs = [("tot", ov[pred[0]].values + ov[pred[1]].values, ov.s.values)] if structure == "tot" else \
                [("D", ov[pred[0]].values, ov.sD.values), ("N", ov[pred[1]].values, ov.sN.values)]
        for name, x, y in pairs:
            sx = np.zeros((ny, self.ng)); sy = np.zeros((ny, self.ng))
            np.add.at(sx, (yi, self.g), x); np.add.at(sy, (yi, self.g), y)
            self.parts[name] = (x, sx, sy)

    def factors(self, exclude_year=None, years=None):
        out = {}
        for name, (x, sx, sy) in self.parts.items():
            keep = np.ones(len(OV_YEARS), bool)
            if years is not None:
                keep &= np.isin(OV_YEARS, years)
            if exclude_year is not None:
                keep &= np.array(OV_YEARS) != exclude_year
            X = sx[keep].sum(0).reshape(self.nz, self.ns); Y = sy[keep].sum(0).reshape(self.nz, self.ns)
            k0 = Y.sum(0) / np.maximum(X.sum(0), 1e-9)                # regional factor per season
            k00 = Y.sum() / max(X.sum(), 1e-9)
            k0 = np.where(X.sum(0) > 50, k0, k00)
            out[name] = ((Y + ALPHA * k0[None, :]) / (X + ALPHA)).reshape(-1)
            out[name + "_global"] = k0
            out[name + "_all"] = k00
        return out

    def predict(self, f, rows=slice(None)):
        h = 0
        for name, (x, _, _) in self.parts.items():
            h = h + f[name][self.g[rows]] * x[rows]
        return h


def evaluate(model, years=None):
    """Leave-one-year-out: hide VIIRS for a year, rebuild it from MODIS, compare."""
    rec = []
    for y in OV_YEARS:
        rows = np.where(ov.year.values == y)[0]
        f = model.factors(exclude_year=y, years=years)
        p = model.predict(f, rows)
        t = pd.DataFrame({"l1": ov.l1.values[rows], "month": ov.month.values[rows], "pred": p, "act": ov.s.values[rows],
                          "mod": (ov.m.values[rows])})
        t = t[t.l1 >= 0].groupby(["l1", "month"], as_index=False).sum()
        t["year"] = y
        rec.append(t)
    return pd.concat(rec, ignore_index=True)


def score(t, minn=20):
    u = t[(t.act >= minn) | (t.pred >= minn)]
    le = np.log((u.pred + 1) / (u.act + 1))
    ann = t.groupby(["l1", "year"], as_index=False)[["pred", "act"]].sum()
    ann = ann[ann.act >= 50]
    la = np.log((ann.pred + 1) / (ann.act + 1))
    return dict(n=int(len(u)), rmsle=float(np.sqrt((le ** 2).mean())), mdape=float(np.median(np.abs(u.pred - u.act) / np.maximum(u.act, 1)) * 100),
                bias_pct=float((u.pred.sum() / u.act.sum() - 1) * 100),
                r2_log=float(np.corrcoef(np.log1p(u.pred), np.log1p(u.act))[0, 1] ** 2),
                annual_mdape=float(np.median(np.abs(ann.pred - ann.act) / ann.act) * 100),
                annual_rmsle=float(np.sqrt((la ** 2).mean())), n_annual=int(len(ann)))


print("\nModel selection (leave-one-year-out, province x month):")
cands = []
for zdeg in (None, 2.0, 1.0, 0.5, 0.3, 0.2):
    for season in ("none", "s4", "m12", "int4", "s4int"):
        for structure in ("tot", "dn"):
            m = Model(zdeg, season, structure)
            t = evaluate(m)
            sc = score(t)
            sc.update(zone="regional" if zdeg is None else "%.1f deg" % zdeg, season=season, structure=structure)
            cands.append((sc, m, t))
            print("  zone %-8s season %-4s %-3s  rmsle %.3f  mdape %5.1f%%  bias %+5.1f%%  annual mdape %5.1f%%" % (
                sc["zone"], season, structure, sc["rmsle"], sc["mdape"], sc["bias_pct"], sc["annual_mdape"]))
best_sc, best, best_t = min(cands, key=lambda c: c[0]["rmsle"])
print("-> chosen: zone %s, season %s, structure %s" % (best_sc["zone"], best_sc["season"], best_sc["structure"]))

# naive baselines for the same test: (a) raw MODIS as if it were VIIRS, (b) one global multiplier
naive_t = best_t.copy(); naive_t["pred"] = naive_t["mod"]
sc_naive = score(naive_t)

# ------------------------------------------------------------------ drift diagnostic + NOAA-20 floor
drift = []
for y in OV_YEARS:
    q = ov[ov.year == y]
    drift.append(dict(year=y, ratio=float(q.s.sum() / max(q.m.sum(), 1)), ratio_day=float(q.sD.sum() / max(q.mD.sum(), 1)),
                      ratio_night=float(q.sN.sum() / max(q.mN.sum(), 1)),
                      terra_share=float(q.t.sum() / max(q.m.sum(), 1)), viirs=int(q.s.sum()), modis=int(q.m.sum())))
floor = None
cJ = {"D": 1.0, "N": 1.0}
if "J" in raw:
    sj = raw["S"][["day", "cell", "nD", "nN"]].merge(raw["J"][["day", "cell", "nD", "nN"]], on=["day", "cell"], how="outer",
                                                      suffixes=("_s", "_j")).fillna(0)
    J0 = int(raw["J"].day.min())
    sj = sj[(sj.day >= J0) & (sj.day <= min(S_END, int(raw["J"].day.max())))]
    cJ = {"D": float(sj.nD_s.sum() / max(sj.nD_j.sum(), 1)), "N": float(sj.nN_s.sum() / max(sj.nN_j.sum(), 1))}
    sj = sj.assign(l1=l1_of_cell[sj.cell.values], year=year_of[sj.day.values], month=month_of[sj.day.values],
                   act=sj.nD_s + sj.nN_s, pred=cJ["D"] * sj.nD_j + cJ["N"] * sj.nN_j)
    tj = sj[sj.l1 >= 0].groupby(["l1", "year", "month"], as_index=False)[["pred", "act"]].sum()
    floor = score(tj)
    floor.update(cJ_day=cJ["D"], cJ_night=cJ["N"], first=d2s(J0))
    print("NOAA-20 vs Suomi NPP (same instrument, 50 min apart): rmsle %.3f mdape %.1f%%  factors %.3f / %.3f" % (
        floor["rmsle"], floor["mdape"], cJ["D"], cJ["N"]))
print("Naive (MODIS counts used as-is): rmsle %.3f mdape %.1f%% bias %+.1f%%" % (sc_naive["rmsle"], sc_naive["mdape"], sc_naive["bias_pct"]))

# ------------------------------------------------------------------ build the harmonised record
f_all = best.factors()
# Terra-only factors (same structure/zone/season, Terra as predictor)
best_T = Model(best.zdeg, best.season, best.structure, pred=("tD", "tN"))
f_T = best_T.factors()


def apply_model(model, f, df, colD, colN):
    c = df.cell.values
    cx, cy = c % NX, c // NX
    if model.zdeg is None:
        z = np.zeros(len(df), np.int32)
    else:
        n = int(round(model.zdeg / 0.1)); z = (cy // n) * ((NX + n - 1) // n) + cx // n
    s, _ = season_of(model.season, c, df.day.values, df[colD].values + df[colN].values, model.scale)
    g = z * model.ns + s
    if model.structure == "tot":
        return f["tot"][g] * (df[colD].values + df[colN].values)
    return f["D"][g] * df[colD].values + f["N"][g] * df[colN].values


pre = M[M.day < S0].copy()
both = pre.day >= A0
pre["h"] = 0.0
pre.loc[both, "h"] = apply_model(best, f_all, pre[both], "mD", "mN")
pre.loc[~both, "h"] = apply_model(best_T, f_T, pre[~both], "T_nD", "T_nN")
post = raw["S"][["day", "cell"]].copy(); post["h"] = (raw["S"].nD + raw["S"].nN).values
H = pd.concat([pre[["day", "cell", "h"]], post], ignore_index=True)
# after Suomi NPP ends, continue with NOAA-20 scaled to Suomi NPP units
if "J" in raw and int(raw["J"].day.max()) > S_END:
    tail = raw["J"][raw["J"].day > S_END]
    H = pd.concat([H, pd.DataFrame({"day": tail.day.values, "cell": tail.cell.values,
                                    "h": cJ["D"] * tail.nD.values + cJ["N"] * tail.nN.values})], ignore_index=True)
Mtot = M[["day", "cell"]].copy(); Mtot["v"] = (M.mD + M.mN).values

# region membership matrix (cells x regions)
rows_, cols_ = [], []
for r in regions:
    c = rc[str(r["id"])]
    rows_.append(c); cols_.append(np.full(len(c), r["id"]))
A = sparse.csr_matrix((np.ones(sum(map(len, rows_)), np.float32), (np.concatenate(rows_), np.concatenate(cols_))), shape=(NCELL, NR))


def region_daily(df, col):
    X = sparse.csr_matrix((df[col].values.astype(np.float32), (df.day.values, df.cell.values)), shape=(D, NCELL))
    return np.asarray((X @ A).todense(), dtype=np.float32).T        # (NR, D)


RH = region_daily(H, "h")
RM = region_daily(Mtot, "v")
print("\nregion-day matrices:", RH.shape, " harmonised total %.0f" % RH[[r["id"] for r in regions if r["name"] == "Whole region"][0]].sum())

# ------------------------------------------------------------------ uncertainty per region (from LOYO annual errors)
f_cache = {y: best.factors(exclude_year=y) for y in OV_YEARS}
pred_loyo = np.zeros(len(ov))
for y in OV_YEARS:
    rows = np.where(ov.year.values == y)[0]
    pred_loyo[rows] = best.predict(f_cache[y], rows)
ny = len(OV_YEARS)
yi = ov.year.values - OV_YEARS[0]
Xp = sparse.csr_matrix((pred_loyo, (yi, ov.cell.values)), shape=(ny, NCELL)) @ A
Xa = sparse.csr_matrix((ov.s.values, (yi, ov.cell.values)), shape=(ny, NCELL)) @ A
Xp, Xa = np.asarray(Xp.todense()), np.asarray(Xa.todense())       # (years, regions)
sig = np.full(NR, np.nan)
lr_all, n_all = [], []
for r in regions:
    i = r["id"]
    ok = Xa[:, i] >= 30
    if ok.sum() >= 6:
        lr = np.log((Xp[ok, i] + 1) / (Xa[ok, i] + 1))
        sig[i] = float(np.sqrt((lr ** 2).mean()))
    okk = Xa[:, i] >= 5
    lr_all += list(np.log((Xp[okk, i] + 1) / (Xa[okk, i] + 1))); n_all += list(Xa[okk, i])
# pooled fallback as a function of size: sigma^2 = a + b / n
lr_all, n_all = np.array(lr_all), np.array(n_all)
bins = [5, 15, 40, 100, 300, 1000, 1e9]
pool = [(float(np.median(n_all[(n_all >= lo) & (n_all < hi)])), float(np.sqrt((lr_all[(n_all >= lo) & (n_all < hi)] ** 2).mean())))
        for lo, hi in zip(bins[:-1], bins[1:]) if ((n_all >= lo) & (n_all < hi)).sum() > 10]
print("pooled sigma_log by annual size:", [(round(a), round(b, 2)) for a, b in pool])


def pooled_sigma(n):
    if not pool:
        return 0.5
    xs, ys = zip(*pool)
    return float(np.interp(np.log(max(n, 1)), np.log(xs), ys))


for r in regions:
    i = r["id"]
    if np.isnan(sig[i]):
        yrs = max(1, (S_END - S0) / 365.25)
        sig[i] = pooled_sigma(RH[i, S0:].sum() / yrs)

# Terra-only era (before Aqua): how much worse is the rebuild when only Terra is available?
sc_T = score(evaluate(best_T))
terra_inflation = float(max(1.0, sc_T["annual_rmsle"] / max(best_sc["annual_rmsle"], 1e-9)))
print("Terra-only rebuild: monthly mdape %.1f%% annual mdape %.1f%%  -> uncertainty x%.2f for 2001" % (sc_T["mdape"], sc_T["annual_mdape"], terra_inflation))
loyo = {str(r["id"]): dict(pred=[round(float(v), 1) for v in Xp[:, r["id"]]], act=[int(v) for v in Xa[:, r["id"]]]) for r in regions if r["level"] <= 1}

# ------------------------------------------------------------------ headline statistics
# a year counts as complete when the record reaches at least 20 December of it
full_years = [y for y in range(2001, int(year_of[-1]) + 1) if (dt.date(y, 12, 20) - EPOCH).days <= D - 1]
ycol = {y: np.where(year_of == y)[0] for y in full_years}


def annual(mat):
    return np.stack([mat[:, ycol[y]].sum(1) for y in full_years], 1)       # (NR, nyears)


AH, AM = annual(RH), annual(RM)
naive = np.where(np.array(full_years)[None, :] >= 2012, AH, AM)             # switch sensors in 2012, no correction


def trend(y, x):
    ok = y > 0
    if ok.sum() < 8:
        return None
    sl = stats.theilslopes(np.log(y[ok]), np.array(x)[ok])
    kt = stats.kendalltau(np.array(x)[ok], y[ok])
    return dict(pct_per_year=float((np.exp(sl[0]) - 1) * 100), lo=float((np.exp(sl[2]) - 1) * 100), hi=float((np.exp(sl[3]) - 1) * 100), p=float(kt.pvalue))


head = {}
for r in regions:
    if r["level"] <= 1:
        i = r["id"]
        pre_m = [k for k, y in enumerate(full_years) if 2003 <= y <= 2011]
        post_m = [k for k, y in enumerate(full_years) if y >= 2012]
        head[r["name"]] = dict(
            trend_harmonised=trend(AH[i], full_years), trend_naive=trend(naive[i], full_years),
            step_naive=float(naive[i, post_m].mean() / max(naive[i, pre_m].mean(), 1e-9)),
            step_harmonised=float(AH[i, post_m].mean() / max(AH[i, pre_m].mean(), 1e-9)),
            top_years=[int(full_years[k]) for k in np.argsort(-AH[i])[:5]], sigma_log=float(sig[i]))
w = head["Whole region"]
print("\nWHOLE REGION  naive step at 2012: x%.2f   harmonised: x%.2f" % (w["step_naive"], w["step_harmonised"]))
print("  trend naive %s\n  trend harmonised %s\n  top years %s" % (w["trend_naive"], w["trend_harmonised"], w["top_years"]))

# ------------------------------------------------------------------ write web data
def enc(v):
    return np.minimum(np.round(np.sqrt(v) * 64), 65535).astype(np.uint16)


def write_series(path, ids):
    index = {}
    with open(path, "wb") as f:
        for i in ids:
            ent = []
            for mat in (RH, RM):
                nz = np.nonzero(mat[i] > 0.02)[0].astype(np.uint16)
                ent.append([f.tell(), int(len(nz))])
                f.write(nz.tobytes()); f.write(enc(mat[i][nz]).tobytes())
            index[str(i)] = ent
    return index


idx_main = write_series(os.path.join(OUT, "series_main.bin"), [r["id"] for r in regions if r["level"] <= 1])
idx_l2 = write_series(os.path.join(OUT, "series_l2.bin"), [r["id"] for r in regions if r["level"] == 2])

# daily 0.2 degree grids for user-drawn areas: offsets per day (uint32) then (cell uint16, value uint16) pairs
NX2, NY2 = 125, 68
def write_daily(path, df, col):
    c = df.cell.values
    c2 = (c // NX // 2) * NX2 + (c % NX) // 2
    t = pd.DataFrame({"day": df.day.values, "c2": c2, "v": df[col].values}).groupby(["day", "c2"], as_index=False).v.sum()
    t = t[t.v > 0.02].sort_values(["day", "c2"])
    off = np.zeros(D + 1, np.uint32)
    np.add.at(off, t.day.values + 1, 1)
    off = np.cumsum(off).astype(np.uint32)
    rec = np.empty((len(t), 2), np.uint16)
    rec[:, 0] = t.c2.values; rec[:, 1] = enc(t.v.values)
    with open(path, "wb") as f:
        f.write(off.tobytes()); f.write(rec.tobytes())
    return len(t)
n_h = write_daily(os.path.join(OUT, "grid_h.bin"), H, "h")
n_m = write_daily(os.path.join(OUT, "grid_m.bin"), Mtot, "v")
land2 = np.zeros(NX2 * NY2, np.uint8)
lc = np.where(land)[0]; land2[(lc // NX // 2) * NX2 + (lc % NX) // 2] = 1

# annual 0.1 degree grid of harmonised activity, log-scaled to uint8
all_years = list(range(int(year_of[0]) + (0 if month_of[0] == 1 else 1), int(year_of[-1]) + 1))
all_years = [y for y in range(2001, int(year_of[-1]) + 1)]
G = np.zeros((len(all_years), NCELL), np.float32)
hy = year_of[H.day.values]
for k, y in enumerate(all_years):
    q = H[hy == y]
    np.add.at(G[k], q.cell.values, q.h.values)
GMAX = float(np.percentile(G[G > 0], 99.9)) if (G > 0).any() else 1.0
G8 = np.where(G > 0, np.clip(1 + np.round(254 * np.log1p(G) / np.log1p(GMAX)), 1, 255), 0).astype(np.uint8)
G8.tofile(os.path.join(OUT, "grid_annual.bin"))

# evidence payload
sc_pts = best_t[(best_t.act >= 5) | (best_t.pred >= 5)]
if len(sc_pts) > 2500:
    sc_pts = sc_pts.sample(2500, random_state=1)
kmap = None
if best.zdeg is not None:
    n = int(round(best.zdeg / 0.1)); nzx = (NX + n - 1) // n; nzy = (NY + n - 1) // n
    name = "tot" if best.structure == "tot" else "D"
    kk = f_all[name].reshape(best.nz, best.ns)
    # only zones with enough MODIS evidence
    x_sum = best.parts[name][1].sum(0).reshape(best.nz, best.ns)
    kmap = dict(n=n, nzx=nzx, nzy=nzy, k=[[round(float(v), 2) for v in row] for row in kk],
                support=[[int(v) for v in row] for row in x_sum])
manifest = json.load(open(os.path.join(DATA, "manifest.json"))) if os.path.exists(os.path.join(DATA, "manifest.json")) else {}
for r in regions:
    r["sigma_log"] = round(float(sig[r["id"]]), 3)
meta = dict(
    built_utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    synthetic=bool(os.environ.get("SYNTHETIC")),
    grid=dict(west=94.5, east=119.5, south=-6.0, north=7.5, res=0.1, nx=NX, ny=NY),
    epoch=str(EPOCH), ndays=D, last_date=d2s(D - 1), s0=S0, a0=A0, t0=int(raw["T"].day.min()), s_end=S_END,
    s0_date=d2s(S0), a0_date=d2s(A0), t0_date=d2s(raw["T"].day.min()),
    full_years=full_years, grid_years=all_years, grid_max=GMAX,
    regions=regions, index_main=idx_main, index_l2=idx_l2,
    sensors={c: dict(first=d2s(v.day.min()), last=d2s(v.day.max()), detections=v.attrs["n_land"]) for c, v in raw.items()},
    model=dict(chosen=best_sc, candidates=[c[0] for c in cands], naive=sc_naive, floor=floor,
               k_regional={k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in f_all.items() if k.endswith("_global") or k.endswith("_all")},
               k_terra_only={k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in f_T.items() if k.endswith("_all")},
               kmap=kmap, alpha=ALPHA, overlap_years=OV_YEARS),
    drift=drift, sigma_pool=pool, loyo=loyo, loyo_years=OV_YEARS, terra_inflation=terra_inflation, terra_only=sc_T, grid2=dict(nx=NX2, ny=NY2, res=0.2),
    recent=dict(sensors=HAS_RECENT, manifest=(json.load(open(os.path.join(RECENT, "manifest_recent.json"))) if RECENT and os.path.exists(os.path.join(RECENT, "manifest_recent.json")) else None)),
    m_end=M_END, j_fill_from=(S_END + 1 if ("J" in raw and int(raw["J"].day.max()) > S_END) else None), cJ=cJ,
    scatter=dict(pred=[round(float(v), 1) for v in sc_pts.pred], act=[int(v) for v in sc_pts.act], year=[int(v) for v in sc_pts.year]),
    headline=head,
    provenance=dict(source=manifest.get("source"), filters=manifest.get("filters"), created=manifest.get("created_utc"),
                    row_stats=manifest.get("row_stats_by_source"), n_files=len(manifest.get("files", []))),
)
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), separators=(",", ":"))
json.dump(dict(headline=head, model=meta["model"]["chosen"], naive=sc_naive, floor=floor, drift=drift), open("pipeline/report.json", "w"), indent=1)
np.save("pipeline/RH.npy", RH); np.save("pipeline/RM.npy", RM)
for fn in sorted(os.listdir(OUT)):
    print("  %-20s %7.2f MB" % (fn, os.path.getsize(os.path.join(OUT, fn)) / 1e6))
