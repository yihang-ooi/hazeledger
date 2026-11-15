#!/usr/bin/env python3
# SI disclosure: the code, the data analysis and the text of this project were produced with Claude (Anthropic), an AI assistant, under the author's direction. No SI-generated images, video or audio are used; every chart is drawn by code from NASA data.
"""
firms_recent.py  -  bring the record up to the present with the NASA FIRMS API
==============================================================================

The yearly country archive used by firms_fetch.py stops at the last complete
year. This script asks the FIRMS area API for everything since 1 Jan 2026
(standard-processing data where NASA has published it, near-real-time after
that) for the same Sumatra / Peninsular Malaysia / Borneo window, and writes
files in the same format.

Run it (needs a free MAP_KEY from https://firms.modaps.eosdis.nasa.gov/api/map_key/):

    python firms_recent.py YOUR_MAP_KEY

Output: ./hazeledger_recent/   (send the whole folder, zipped, to Claude)
Python 3.8+ standard library only. Responses are cached in ./firms_cache/.
"""
import csv
import datetime as dt
import gzip
import io
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.request

API = os.environ.get("FIRMS_API", "https://firms.modaps.eosdis.nasa.gov/api")
START = dt.date(2026, 1, 1)
TODAY = dt.datetime.now(dt.timezone.utc).date()
CACHE_DIR = "firms_cache"
OUT_DIR = "hazeledger_recent"
WEST, EAST, SOUTH, NORTH = 94.5, 119.5, -6.0, 7.5
AREA = "%g,%g,%g,%g" % (WEST, SOUTH, EAST, NORTH)
RES = 0.1
NX, NY = 250, 135
NCELL = NX * NY
EPOCH = dt.date(2000, 1, 1)
CHUNK = [5]   # days per request; shrinks automatically if NASA refuses the size

# family -> (standard-processing source, near-real-time source)
FAMILIES = {
    "MODIS": ("MODIS_SP", "MODIS_NRT"),
    "SNPP": ("VIIRS_SNPP_SP", "VIIRS_SNPP_NRT"),
    "NOAA20": ("VIIRS_NOAA20_SP", "VIIRS_NOAA20_NRT"),
    "NOAA21": ("VIIRS_NOAA21_SP", "VIIRS_NOAA21_NRT"),
}
SAT_CODE = {"TERRA": "T", "AQUA": "A", "N": "S", "NPP": "S", "SNPP": "S",
            "N20": "J", "NOAA-20": "J", "N21": "K", "NOAA-21": "K"}
DEFAULT_CODE = {"SNPP": "S", "NOAA20": "J", "NOAA21": "K"}
CODE_NAME = {"T": "MODIS Terra", "A": "MODIS Aqua", "S": "VIIRS Suomi NPP", "J": "VIIRS NOAA-20", "K": "VIIRS NOAA-21"}
UA = {"User-Agent": "Mozilla/5.0 (SpaceApps2026 student project; FIRMS area API)"}


def ssl_context():
    try:
        import certifi  # type: ignore
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


CTX = ssl_context()


class Fatal(Exception):
    pass


class BadRequest(str):
    """Marker for an HTTP 400/404 answer; carries NASA's message."""
    def __new__(cls, code, body):
        o = str.__new__(cls, "")
        o.code, o.body = code, body
        return o


DEAD = set()   # sources the API refused


def http_text(url, tries=8):
    last = ""
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=120, context=CTX) as r:
                txt = r.read().decode("utf-8", "replace")
            low = txt[:300].lower()
            if "invalid map_key" in low or "invalid map key" in low:
                raise Fatal("NASA says the MAP_KEY is invalid. Check for a typo or request a new key.")
            if "exceed" in low and "transaction" in low:
                print("    NASA rate limit reached, waiting 90 s ...")
                sys.stdout.flush()
                time.sleep(90)
                continue
            return txt
        except Fatal:
            raise
        except urllib.error.HTTPError as e:
            last = "HTTP %s" % e.code
            if e.code in (429, 403):
                print("    NASA rate limit (HTTP %s), waiting 90 s ..." % e.code)
                sys.stdout.flush()
                time.sleep(90)
                continue
            if e.code in (400, 404):
                try:
                    body = e.read().decode("utf-8", "replace").strip()[:160]
                except Exception:
                    body = ""
                return BadRequest(e.code, body)
        except Exception as e:
            last = "%s: %s" % (type(e).__name__, e)
        time.sleep(4 * (attempt + 1))
    print("    giving up on one request (%s)" % last)
    return None


def availability(key):
    """source id -> (min_date, max_date); empty dict if NASA's answer cannot be parsed."""
    out = {}
    txt = http_text("%s/data_availability/csv/%s/ALL" % (API, key), tries=3) or ""
    try:
        rows = list(csv.reader(io.StringIO(txt)))
        for row in rows[1:]:
            ds = [c.strip() for c in row]
            dates = [c for c in ds[1:] if len(c) == 10 and c[4] == "-" and c[7] == "-"]
            if len(dates) >= 2:
                out[ds[0]] = (dt.date.fromisoformat(dates[0]), dt.date.fromisoformat(dates[1]))
    except Exception:
        pass
    return out


def fetch_chunk(key, source, start, ndays):
    """Returns list of dict rows (possibly empty), or None on failure. Cached on disk when complete."""
    if source in DEAD:
        return []
    name = os.path.join(CACHE_DIR, "api_%s_%s_%d.csv" % (source, start.isoformat(), ndays))
    final = start + dt.timedelta(days=ndays + 2) < TODAY     # only cache windows that can no longer change
    if final and os.path.exists(name):
        txt = open(name, encoding="utf-8").read()
    else:
        txt = http_text("%s/area/csv/%s/%s/%s/%d/%s" % (API, key, source, AREA, ndays, start.isoformat()))
        if txt is None:
            return None
        if isinstance(txt, BadRequest):
            return txt
        if txt.strip() and not txt.lstrip().lower().startswith("latitude"):
            if source not in DEAD:
                print("    NASA answered %r for %s; skipping this source" % (txt.strip()[:80], source))
            DEAD.add(source)
            return []
        if final:
            with open(name, "w", encoding="utf-8") as f:
                f.write(txt)
        time.sleep(0.25)
    if not txt.strip():
        return []
    return list(csv.DictReader(io.StringIO(txt)))


def main():
    if len(sys.argv) < 2 or len(sys.argv[1].strip()) < 20:
        print(__doc__)
        print("ERROR: pass your MAP_KEY, for example:  python firms_recent.py 0123456789abcdef0123456789abcdef")
        sys.exit(1)
    key = sys.argv[1].strip()
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    avail = availability(key)
    print("NASA data availability:" if avail else "Could not read NASA's availability table; probing instead.")
    for k in sorted(avail):
        print("   %-18s %s to %s" % (k, avail[k][0], avail[k][1]))

    aggs = {}          # code -> {key: [nD, nN, frpD, frpN]}
    static_cells = {}  # cell -> count of non-vegetation detections (standard data only)
    nrt_rows = []      # (code, day, cell, isNight, frp) kept aside until the static mask is known
    info = {}
    for fam, (sp, nrt) in FAMILIES.items():
        sp_max = avail.get(sp, (None, None))[1] if avail else None
        nrt_rng = avail.get(nrt) if avail else None
        if avail and sp not in avail and nrt not in avail:
            print("\n%s: not offered by the API, skipped" % fam)
            continue
        print("\n%s" % fam)
        st = dict(sp_rows=0, nrt_rows=0, sp_last=None, nrt_first=None, nrt_last=None, failed_windows=0)
        d = START
        while d <= TODAY:
            n = min(CHUNK[0], (TODAY - d).days + 1)
            end = d + dt.timedelta(days=n - 1)
            rows_sp, rows_nrt = [], []
            if sp_max is None or d <= sp_max:
                rows_sp = fetch_chunk(key, sp, d, n)
                if isinstance(rows_sp, BadRequest):
                    if CHUNK[0] > 1:
                        CHUNK[0] = {5: 3, 3: 2, 2: 1}.get(CHUNK[0], 1)
                        print("   NASA refused the request (HTTP %s %s); retrying with %d-day windows" % (rows_sp.code, rows_sp.body, CHUNK[0]))
                        continue
                    print("   NASA refused %s for %s (HTTP %s %s)" % (sp, d, rows_sp.code, rows_sp.body))
                    rows_sp = None
                if rows_sp is None:
                    st["failed_windows"] += 1
                    rows_sp = []
            last_sp = max((r["acq_date"] for r in rows_sp), default=None)
            need_nrt = (sp_max is not None and end > sp_max) or (sp_max is None and (last_sp is None or last_sp < end.isoformat()))
            if need_nrt and (nrt_rng is None or end >= nrt_rng[0]):
                rows_nrt = fetch_chunk(key, nrt, d, n)
                if isinstance(rows_nrt, BadRequest):
                    print("   NASA refused %s for %s (HTTP %s %s)" % (nrt, d, rows_nrt.code, rows_nrt.body))
                    rows_nrt = None
                if rows_nrt is None:
                    st["failed_windows"] += 1
                    rows_nrt = []
                cut = sp_max.isoformat() if sp_max is not None else (last_sp or "0000-00-00")
                rows_nrt = [r for r in rows_nrt if r["acq_date"] > cut]
            for is_nrt, rows in ((False, rows_sp), (True, rows_nrt)):
                for r in rows:
                    try:
                        lat, lon = float(r["latitude"]), float(r["longitude"])
                        if not (SOUTH <= lat < NORTH and WEST <= lon < EAST):
                            continue
                        ix = min(NX - 1, int((lon - WEST) / RES)); iy = min(NY - 1, int((lat - SOUTH) / RES))
                        cell = iy * NX + ix
                        typ = (r.get("type") or "").strip()
                        if typ not in ("", "0"):
                            static_cells[cell] = static_cells.get(cell, 0) + 1
                            continue
                        c = r["confidence"].strip().lower()
                        if fam == "MODIS":
                            if float(c) < 30:
                                continue
                        elif c in ("l", "low"):
                            continue
                        code = SAT_CODE.get((r.get("satellite") or "").strip().upper(), DEFAULT_CODE.get(fam))
                        if code is None:
                            continue
                        day = (dt.date.fromisoformat(r["acq_date"].strip()) - EPOCH).days
                        night = (r.get("daynight") or "D").strip().upper().startswith("N")
                        frp = float(r.get("frp") or 0)
                    except (ValueError, KeyError):
                        continue
                    if is_nrt:
                        nrt_rows.append((code, day, cell, night, frp))
                        st["nrt_rows"] += 1
                        st["nrt_first"] = min(st["nrt_first"] or r["acq_date"], r["acq_date"])
                        st["nrt_last"] = max(st["nrt_last"] or r["acq_date"], r["acq_date"])
                    else:
                        add(aggs, code, day, cell, night, frp)
                        st["sp_rows"] += 1
                        st["sp_last"] = max(st["sp_last"] or r["acq_date"], r["acq_date"])
            print("   %s to %s   standard %6d   near-real-time %6d" % (d, end, len(rows_sp), len(rows_nrt)))
            sys.stdout.flush()
            d = end + dt.timedelta(days=1)
        info[fam] = st

    # near-real-time rows carry no type flag: drop those sitting on cells that standard data marks as static heat sources
    mask = {c for c, n in static_cells.items() if n >= 3}
    dropped = 0
    for code, day, cell, night, frp in nrt_rows:
        if cell in mask:
            dropped += 1
            continue
        add(aggs, code, day, cell, night, frp)

    for code, a in sorted(aggs.items()):
        with gzip.open(os.path.join(OUT_DIR, "recent_%s.csv.gz" % code), "wt", compresslevel=6) as w:
            w.write("day,cell,nD,nN,frpD,frpN\n")
            for k in sorted(a):
                nD, nN, fD, fN = a[k]
                w.write("%d,%d,%d,%d,%.1f,%.1f\n" % (k // NCELL, k % NCELL, nD, nN, fD, fN))
    manifest = dict(created_utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    source="NASA FIRMS area API, %s/area/" % API, start=str(START), end=str(TODAY), area=AREA,
                    families=info, availability={k: [str(v[0]), str(v[1])] for k, v in avail.items()},
                    static_cells_masked=len(mask), nrt_rows_dropped_as_static=dropped,
                    filters={"type": "standard data: type 0 only; near-real-time: static-source cells masked",
                             "confidence": "MODIS >= 30; VIIRS nominal + high"},
                    sensors={c: {"name": CODE_NAME[c], "first_date": str(EPOCH + dt.timedelta(days=min(a) // NCELL)),
                                 "last_date": str(EPOCH + dt.timedelta(days=max(a) // NCELL)),
                                 "detections": sum(v[0] + v[1] for v in a.values())} for c, a in aggs.items() if a})
    with open(os.path.join(OUT_DIR, "manifest_recent.json"), "w") as f:
        json.dump(manifest, f, indent=1)

    print("\n" + "=" * 70)
    print("DONE in %.0f s" % (time.time() - t0))
    for c, s in sorted(manifest["sensors"].items()):
        print("  %-16s %9d detections   %s to %s" % (s["name"], s["detections"], s["first_date"], s["last_date"]))
    bad = sum(v["failed_windows"] for v in info.values())
    if bad:
        print("\n%d request(s) failed. Run the script again; finished windows are cached." % bad)
    for fam, st in info.items():
        if fam != "NOAA21" and st["sp_rows"] == 0:
            print("\n*** WARNING: no standard-processing rows came back for %s. Paste this whole output to Claude." % fam)
    print("\nSend ALL files in the folder %r to Claude:" % OUT_DIR)
    for fn in sorted(os.listdir(OUT_DIR)):
        print("   %-24s %6.2f MB" % (fn, os.path.getsize(os.path.join(OUT_DIR, fn)) / 1e6))
    print("=" * 70)


def add(aggs, code, day, cell, night, frp):
    a = aggs.setdefault(code, {})
    rec = a.get(day * NCELL + cell)
    if rec is None:
        rec = [0, 0, 0.0, 0.0]
        a[day * NCELL + cell] = rec
    if night:
        rec[1] += 1; rec[3] += frp
    else:
        rec[0] += 1; rec[2] += frp


if __name__ == "__main__":
    try:
        main()
    except Fatal as e:
        print("\nSTOPPED: %s" % e)
        sys.exit(1)
