#!/usr/bin/env python3
"""
firms_fetch.py  -  NASA FIRMS archive downloader + aggregator
=============================================================

What it does
------------
1. Downloads NASA FIRMS "yearly summary by country" active-fire CSV files
   (https://firms.modaps.eosdis.nasa.gov/country/) for Indonesia, Malaysia,
   Brunei and Singapore, for three sensors:
       modis        MODIS Terra + Aqua, 1 km   (MCD14ML, 2000 onward)
       viirs-snpp   VIIRS Suomi NPP, 375 m     (VNP14IMGT, 2012 onward)
       viirs-jpss1  VIIRS NOAA-20, 375 m       (VJ114IMGT, 2018 onward)
   (NOAA-21 is tried too and skipped quietly if NASA has not published it.)
2. Keeps only vegetation-fire detections (type 0) at nominal or high
   confidence inside the Sumatra / Peninsular Malaysia / Borneo window.
3. Counts detections per day, per 0.1 degree grid cell, per satellite,
   split into day and night overpasses, and sums fire radiative power.
4. Writes a handful of small files into  ./hazeledger_data/

Run it
------
    python firms_fetch.py

No extra packages needed (Python 3.8+ standard library only).
Raw NASA files are cached in ./firms_cache/ so a re-run does not download again.
Expect roughly 1.5 to 2.5 GB of downloads in total.
"""

import csv
import datetime as dt
import gzip
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
BASE = os.environ.get("FIRMS_BASE", "https://firms.modaps.eosdis.nasa.gov/data/country")
CACHE_DIR = "firms_cache"
OUT_DIR = "hazeledger_data"
THIS_YEAR = dt.date.today().year

# (folder name on the NASA server, first year to try)
SOURCES = [
    ("modis", 2000),
    ("viirs-snpp", 2012),
    ("viirs-jpss1", 2018),
    ("viirs-jpss2", 2023),
]

# Country -> possible spellings used in NASA file names (first hit wins)
COUNTRIES = {
    "Indonesia": ["Indonesia"],
    "Malaysia": ["Malaysia"],
    "Brunei": ["Brunei", "Brunei_Darussalam", "Brunei%20Darussalam"],
    "Singapore": ["Singapore"],
}

# Study window and grid (0.1 degree cells)
WEST, EAST, SOUTH, NORTH = 94.5, 119.5, -6.0, 7.5
RES = 0.1
NX = int(round((EAST - WEST) / RES))   # 250
NY = int(round((NORTH - SOUTH) / RES))  # 135
NCELL = NX * NY
EPOCH = dt.date(2000, 1, 1)

# Satellite name in the CSV -> one-letter code used in output file names
SAT_CODE = {"TERRA": "T", "AQUA": "A", "N": "S", "NPP": "S", "SNPP": "S",
            "N20": "J", "NOAA-20": "J", "N21": "K", "NOAA-21": "K"}
FOLDER_DEFAULT_CODE = {"viirs-snpp": "S", "viirs-jpss1": "J", "viirs-jpss2": "K"}
CODE_NAME = {"T": "MODIS Terra", "A": "MODIS Aqua", "S": "VIIRS Suomi NPP",
             "J": "VIIRS NOAA-20", "K": "VIIRS NOAA-21"}

UA = {"User-Agent": "Mozilla/5.0 (SpaceApps2026 student project; FIRMS country archive)"}


def ssl_context():
    """Use certifi's certificate bundle when available (fixes macOS python.org installs)."""
    try:
        import certifi  # type: ignore
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


CTX = ssl_context()


# ----------------------------------------------------------------------------
# Download
# ----------------------------------------------------------------------------
def fetch(url, dest, tries=4):
    """Download url to dest. Returns (status, bytes). status: ok | cached | missing | error:<why>"""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached", os.path.getsize(dest)
    part = dest + ".part"
    last = ""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=90, context=CTX) as r, open(part, "wb") as f:
                n = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    n += len(chunk)
            with open(part, "rb") as f:
                head = f.read(200).lstrip().lower()
            if not head.startswith(b"latitude"):
                os.remove(part)
                return "missing", 0  # server answered with a web page, not a CSV
            os.replace(part, dest)
            return "ok", n
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 410):
                return "missing", 0
            last = "HTTP %s" % e.code
        except ssl.SSLError as e:
            return "error:SSL certificate problem (%s)" % e, 0
        except Exception as e:  # timeouts, resets
            last = "%s: %s" % (type(e).__name__, e)
            if "CERTIFICATE_VERIFY_FAILED" in last:
                return "error:SSL certificate problem", 0
        time.sleep(3 * (attempt + 1))
    if os.path.exists(part):
        os.remove(part)
    return "error:" + last, 0


def job_list():
    jobs = []
    for folder, first in SOURCES:
        for year in range(first, THIS_YEAR + 1):
            for country, spellings in COUNTRIES.items():
                jobs.append((folder, year, country, spellings))
    return jobs


def download_one(job):
    folder, year, country, spellings = job
    for name in spellings:
        fname = "%s_%d_%s.csv" % (folder, year, name)
        url = "%s/%s/%d/%s" % (BASE, folder, year, fname)
        dest = os.path.join(CACHE_DIR, fname)
        status, nbytes = fetch(url, dest)
        if status in ("ok", "cached") or status.startswith("error"):
            return dict(folder=folder, year=year, country=country, url=url,
                        file=dest, status=status, bytes=nbytes)
    return dict(folder=folder, year=year, country=country, url=url,
                file=None, status="missing", bytes=0)


# ----------------------------------------------------------------------------
# Parse + aggregate
# ----------------------------------------------------------------------------
_day_cache = {}


def day_index(s):
    d = _day_cache.get(s)
    if d is None:
        d = (dt.date(int(s[0:4]), int(s[5:7]), int(s[8:10])) - EPOCH).days
        _day_cache[s] = d
    return d


def parse_file(path, folder, aggs, stats):
    """Stream one FIRMS CSV into aggs[code][day*NCELL+cell] = [nD, nN, frpD, frpN]."""
    is_modis = folder == "modis"
    with open(path, "r", newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f)
        try:
            header = [h.strip().lower() for h in next(rd)]
        except StopIteration:
            return
        col = {h: i for i, h in enumerate(header)}
        need = ["latitude", "longitude", "acq_date", "confidence", "frp", "daynight"]
        for k in need:
            if k not in col:
                raise ValueError("column %r not found in %s (header: %s)" % (k, path, header))
        i_lat, i_lon, i_date = col["latitude"], col["longitude"], col["acq_date"]
        i_conf, i_frp, i_dn = col["confidence"], col["frp"], col["daynight"]
        i_sat = col.get("satellite")
        i_type = col.get("type")
        default_code = FOLDER_DEFAULT_CODE.get(folder)
        for row in rd:
            stats["rows_read"] += 1
            try:
                lat = float(row[i_lat])
                lon = float(row[i_lon])
                if not (SOUTH <= lat < NORTH and WEST <= lon < EAST):
                    stats["outside_window"] += 1
                    continue
                if i_type is not None:
                    t = row[i_type].strip()
                    if t not in ("", "0"):
                        stats["dropped_type_%s" % t] += 1
                        continue
                c = row[i_conf].strip().lower()
                if is_modis:
                    if float(c) < 30:
                        stats["dropped_lowconf"] += 1
                        continue
                elif c in ("l", "low"):
                    stats["dropped_lowconf"] += 1
                    continue
                sat = row[i_sat].strip().upper() if i_sat is not None else ""
                code = SAT_CODE.get(sat, default_code)
                if code is None:
                    stats["unknown_satellite"] += 1
                    continue
                day = day_index(row[i_date].strip())
                ix = int((lon - WEST) / RES)
                iy = int((lat - SOUTH) / RES)
                if ix >= NX:
                    ix = NX - 1
                if iy >= NY:
                    iy = NY - 1
                key = day * NCELL + iy * NX + ix
                try:
                    frp = float(row[i_frp])
                except ValueError:
                    frp = 0.0
                a = aggs.setdefault(code, {})
                rec = a.get(key)
                if rec is None:
                    rec = [0, 0, 0.0, 0.0]
                    a[key] = rec
                if row[i_dn].strip().upper().startswith("N"):
                    rec[1] += 1
                    rec[3] += frp
                else:
                    rec[0] += 1
                    rec[2] += frp
                stats["kept"] += 1
            except (ValueError, IndexError):
                stats["bad_rows"] += 1


class Counter(dict):
    def __missing__(self, k):
        return 0


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    jobs = job_list()
    print("NASA FIRMS country archive -> %d candidate files (missing years are skipped)" % len(jobs))
    print("Cache folder: %s    Output folder: %s\n" % (os.path.abspath(CACHE_DIR), os.path.abspath(OUT_DIR)))

    # ---- 1. download ---------------------------------------------------------
    results = []
    t0 = time.time()
    done_bytes = 0
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = [ex.submit(download_one, j) for j in jobs]
        for k, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            results.append(r)
            done_bytes += r["bytes"]
            if r["status"] != "missing":
                print("  [%3d/%d] %-9s %-11s %d %-10s %7.1f MB   (total %.0f MB, %.0f s)" % (
                    k, len(jobs), r["status"][:9], r["folder"], r["year"], r["country"],
                    r["bytes"] / 1e6, done_bytes / 1e6, time.time() - t0))
                if r["status"].startswith("error"):
                    print("        -> " + r["status"])
            sys.stdout.flush()

    errors = [r for r in results if r["status"].startswith("error")]
    good = [r for r in results if r["status"] in ("ok", "cached")]
    if any("SSL" in r["status"] for r in errors):
        print("\n*** SSL certificate problem on this Python install. Easiest fixes:")
        print("    a) run:  pip install certifi   then run this script again, or")
        print("    b) run the script in Google Colab instead (see instructions).")
    if not good:
        print("\nNo files could be downloaded. Please copy everything above and send it to Claude.")
        sys.exit(1)
    for core in ("modis", "viirs-snpp", "viirs-jpss1"):
        if not any(r["folder"] == core and r["country"] == "Indonesia" for r in good):
            print("\n*** WARNING: no Indonesia files found for %r. Tell Claude: this sensor is missing." % core)

    # ---- 2. parse, one sensor-year at a time ----------------------------------
    print("\nParsing and aggregating ...")
    writers = {}
    year_table = {}     # (code, year) -> detections kept
    stats_by_folder = {}
    date_range = {}     # code -> [min_day, max_day]
    groups = {}
    for r in good:
        groups.setdefault((r["folder"], r["year"]), []).append(r)
    for (folder, year) in sorted(groups):
        aggs = {}
        stats = stats_by_folder.setdefault(folder, Counter())
        for r in sorted(groups[(folder, year)], key=lambda x: x["country"]):
            before = stats["kept"]
            parse_file(r["file"], folder, aggs, stats)
            r["rows_kept"] = stats["kept"] - before
        for code, a in sorted(aggs.items()):
            if code not in writers:
                writers[code] = gzip.open(os.path.join(OUT_DIR, "agg_%s.csv.gz" % code), "wt", compresslevel=6)
                writers[code].write("day,cell,nD,nN,frpD,frpN\n")
            w = writers[code]
            n = 0
            for key in sorted(a):
                nD, nN, fD, fN = a[key]
                w.write("%d,%d,%d,%d,%.1f,%.1f\n" % (key // NCELL, key % NCELL, nD, nN, fD, fN))
                n += nD + nN
            year_table[(code, year)] = n
            lo, hi = min(a) // NCELL, max(a) // NCELL
            dr = date_range.setdefault(code, [lo, hi])
            dr[0], dr[1] = min(dr[0], lo), max(dr[1], hi)
        print("  %-11s %d  ->  %s" % (folder, year, ", ".join(
            "%s %d" % (CODE_NAME[c], year_table[(c, year)]) for c in sorted(aggs)) or "nothing inside window"))
        sys.stdout.flush()
    for w in writers.values():
        w.close()

    # ---- 3. extras: ENSO index (NOAA), best effort ----------------------------
    oni_status = fetch_text("https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
                            os.path.join(OUT_DIR, "oni.ascii.txt"))

    # ---- 4. manifest -----------------------------------------------------------
    manifest = {
        "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": "NASA FIRMS yearly summary by country, https://firms.modaps.eosdis.nasa.gov/country/",
        "base_url": BASE,
        "grid": {"west": WEST, "east": EAST, "south": SOUTH, "north": NORTH, "res_deg": RES,
                 "nx": NX, "ny": NY, "cell_index": "cell = iy*nx + ix, iy from south, ix from west",
                 "day_index": "days since 2000-01-01 (acq_date, UTC)"},
        "filters": {"type": "kept type 0 (presumed vegetation fire) only",
                    "confidence": "MODIS >= 30; VIIRS nominal + high"},
        "columns": "day,cell,nD,nN,frpD,frpN  (day/night detection counts, summed FRP in MW)",
        "sensors": {c: {"name": CODE_NAME[c],
                        "first_date": str(EPOCH + dt.timedelta(days=date_range[c][0])),
                        "last_date": str(EPOCH + dt.timedelta(days=date_range[c][1]))}
                    for c in sorted(date_range)},
        "detections_kept_by_sensor_year": {"%s_%d" % k: v for k, v in sorted(year_table.items())},
        "row_stats_by_source": {k: dict(v) for k, v in stats_by_folder.items()},
        "files": [{k: r.get(k) for k in ("folder", "year", "country", "url", "status", "bytes", "rows_kept")}
                  for r in sorted(results, key=lambda x: (x["folder"], x["year"], x["country"]))
                  if r["status"] != "missing"],
        "missing": ["%s %d %s" % (r["folder"], r["year"], r["country"])
                    for r in sorted(results, key=lambda x: (x["folder"], x["year"], x["country"]))
                    if r["status"] == "missing"],
        "oni": oni_status,
    }
    with open(os.path.join(OUT_DIR, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)

    # ---- 5. summary ------------------------------------------------------------
    print("\n" + "=" * 70)
    print("DONE in %.0f s. Detections kept, by satellite:" % (time.time() - t0))
    for c in sorted(date_range):
        tot = sum(v for (cc, _), v in year_table.items() if cc == c)
        print("  %-16s %10d   %s to %s" % (CODE_NAME[c], tot, manifest["sensors"][c]["first_date"],
                                           manifest["sensors"][c]["last_date"]))
    if errors:
        print("\n%d file(s) FAILED to download (run the script again to retry just those):" % len(errors))
        for r in errors[:15]:
            print("   %s %d %s  %s" % (r["folder"], r["year"], r["country"], r["status"]))
    print("\nSend ALL of these files to Claude:")
    for fn in sorted(os.listdir(OUT_DIR)):
        p = os.path.join(OUT_DIR, fn)
        print("   %-22s %6.1f MB" % (fn, os.path.getsize(p) / 1e6))
    print("=" * 70)


def fetch_text(url, dest):
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=40, context=CTX) as r:
            data = r.read()
        with open(dest, "wb") as f:
            f.write(data)
        return "ok"
    except Exception as e:
        return "skipped (%s)" % type(e).__name__


if __name__ == "__main__":
    main()
