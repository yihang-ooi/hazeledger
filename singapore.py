"""Singapore air quality vs the harmonised fire record. Inputs: NEA Historical 24-hr PSI (data.gov.sg) and
NEA annual PM10 24-hr mean 99th percentile (data.gov.sg dataset d_dc29ecdc06f1f2b42c3dd06dd1ab8e6b)."""
import sys, json, warnings, numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
PSI = sys.argv[1]
PM10_P99 = {2000: 64, 2001: 56, 2002: 72, 2003: 52, 2004: 60, 2005: 61, 2006: 120, 2007: 53, 2008: 49, 2009: 59, 2010: 76, 2011: 55, 2012: 57,
            2013: 215, 2014: 75, 2015: 186, 2016: 61, 2017: 57, 2018: 59, 2019: 90, 2020: 43, 2021: 51, 2022: 51, 2023: 62, 2024: 52, 2025: 51}
RH = np.load("pipeline/RH.npy")
m = json.load(open("web/data/meta.json"))
R = {r["name"]: r["id"] for r in m["regions"] if r["level"] <= 1}
D = RH.shape[1]
dates = pd.date_range("2000-01-01", periods=D)
sum_ = pd.Series(RH[R["Sumatra"]], index=dates)
riau = pd.Series(RH[R["Riau"]], index=dates)
p = pd.read_csv(PSI)
t = pd.to_datetime(p["24hr_psi"], format="%d/%m/%Y %H:%M", errors="coerce")
cols = ["north", "south", "east", "west", "central"]
for c in cols:
    p[c] = pd.to_numeric(p[c], errors="coerce")
p["mx"] = p[cols].max(axis=1); p["t"] = t; p = p.dropna(subset=["t", "mx"])
daily = p.groupby(p.t.dt.normalize()).mx.max()
j = pd.concat([daily.rename("psi"), sum_.rename("fire")], axis=1, sort=True).dropna()
mo = j.resample("MS").agg({"psi": "mean", "fire": "sum"})
mo["n"] = j.psi.resample("MS").size(); mo["bad"] = (j.psi > 100).resample("MS").sum()
mo = mo[mo.n >= 25]
j7 = j.rolling(7).mean().dropna()
years = [y for y in m["full_years"]]
ann = [float(sum_[str(y)].sum()) for y in years]
junes = {y: float(riau[(riau.index.year == y) & (riau.index.month == 6)].sum()) for y in years}
rank13 = 1 + sum(v > junes[2013] for v in junes.values())
out = dict(
    months=[d.strftime("%Y-%m") for d in mo.index], fire=[round(float(v)) for v in mo.fire], psi=[round(float(v), 1) for v in mo.psi], bad=[int(v) for v in mo.bad],
    years=years, fire_annual=[round(v) for v in ann], pm10=[PM10_P99[y] for y in years],
    stats=dict(monthly_spearman=float(stats.spearmanr(mo.fire, mo.psi)[0]), n_months=int(len(mo)), weekly_spearman=float(stats.spearmanr(j7.fire, j7.psi)[0]),
               n_days=int(len(j)), annual_spearman=float(stats.spearmanr(ann, [PM10_P99[y] for y in years])[0]),
               psi_first=str(j.index.min().date()), psi_last=str(j.index.max().date()),
               riau_june_2013=round(junes[2013]), riau_june_2013_rank=int(rank13), riau_june_median=round(float(np.median(list(junes.values())))),
               sumatra_2013_rank=int(1 + sum(a > ann[years.index(2013)] for a in ann)),
               worst_psi_day=str(daily.idxmax().date()), worst_psi=int(daily.max())))
json.dump(out, open("web/data/singapore.json", "w"), separators=(",", ":"))
print(json.dumps(out["stats"], indent=1))
