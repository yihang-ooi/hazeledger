# HazeLedger

**26 years of fire, in one record that survives the satellites.**

A harmonised burning activity calendar for Sumatra, Borneo and the Malay Peninsula, built from NASA FIRMS active fire data, November 2000 to 14 November 2026.
Entry for the NASA Space Apps Challenge 2026 (Singapore), challenge: *Harmonization of MODIS and VIIRS Hot Spots*.

Live app: the GitHub Pages link in this repository's About box.

## The problem

NASA has counted fires with MODIS (1 km pixels) since 2000 and with VIIRS (375 m pixels) since 2012.
In this region VIIRS reports about 4.8 times as many detections as MODIS for the same burning.
Stitch the two archives together and every province appears to burn about three times more after 2012.
That jump is the satellite, not the fires. MODIS retires in 2027, so the seam is about to move again.

## What HazeLedger does

- Joins MODIS (Terra, Aqua) and VIIRS (Suomi NPP, NOAA-20) into one daily record in a single unit: VIIRS-equivalent detections.
- Shows it as a burning calendar by day or by month, for 361 places (island groups, provinces, districts) or any rectangle you draw or type as coordinates.
- Flags unusual months, marks the critical burning windows, and compares the latest month with its historical median.
- Lets you pick a date range and download the daily record, the monthly calendar (CSV) and the calendar image (PNG).
- Loads NASA true-colour imagery for any day you click, so you can see the smoke.
- Shows its own evidence: every year of the 2012 to 2025 overlap was hidden in turn and rebuilt from MODIS alone.

## Results on hidden years

| Approach | Typical monthly error, per province |
|---|---|
| MODIS counts used as they are | 80% |
| HazeLedger | 31% (12% on annual totals, R² 0.83 on a log scale) |
| Two identical VIIRS instruments 50 minutes apart (the floor) | 16% |

Raw stitched counts say Sumatra's average year rose 2.9 times after 2012. The harmonised record says it fell to 0.62 times, in step with MODIS alone (0.60 times), which never changed sensors.

## What the record says about 2026

October 2026 brought 119,994 VIIRS-equivalent detections across the region, 9.4 times the median October and the third highest October in 26 years. For Kalimantan, 2026 is already the third worst year on record.

## NASA data

- **NASA FIRMS active fire archive**, yearly summary by country, 2000 to 2025: https://firms.modaps.eosdis.nasa.gov/country/
  - MODIS Collection 6.1 (MCD14ML), Terra and Aqua
  - VIIRS 375 m Suomi NPP (VNP14IMGT)
  - VIIRS 375 m NOAA-20 (VJ114IMGT)
- **NASA FIRMS area API**, 1 January to 14 November 2026 (standard processing to mid 2026, near-real-time after): https://firms.modaps.eosdis.nasa.gov/api/area/
- **NASA GIBS** imagery (MODIS corrected reflectance and thermal anomalies): https://gibs.earthdata.nasa.gov/

5,305,379 detections were kept after filtering: vegetation fires only (type 0; for near-real-time data, cells with known static heat sources are masked), nominal or high confidence, on land inside 94.5°E to 119.5°E and 6°S to 7.5°N.
Calibration and all hold-out tests use the 2012 to 2025 archive years only.

We acknowledge the use of data from NASA's Fire Information for Resource Management System (FIRMS), part of NASA's Earth Science Data and Information System (ESDIS). NASA does not endorse this project.

## Other open data

- Province and district boundaries: geoBoundaries (CC BY 4.0), https://www.geoboundaries.org/
- Singapore air quality: National Environment Agency via data.gov.sg (Historical 24-hr PSI; PM10 24-hr mean 99th percentile)

## Method in brief

1. Count every detection on a 0.1° grid, per day, per satellite.
2. In the overlap years, learn the ratio of VIIRS to MODIS detections for each 0.3° block and each of four burning intensities, shrunk toward the regional ratio where a block has few fires.
3. Before 20 January 2012, multiply MODIS by the local factor. From then on, use VIIRS itself (Suomi NPP to 2 November 2026, then NOAA-20).
4. Choose the model by leave-one-year-out error (60 variants compared) and turn the hold-out errors into the uncertainty range shown on rebuilt years.

Limits: fires too small or too cool for MODIS on a given day are restored in totals, not in timing. Cloud and thick smoke hide fires from both sensors. A hotspot is an observation of heat, not proof of who lit a fire.

## Files

- `index.html` and the `.bin`, `.json`, `.geojson` files beside it: the app and its data.
- `firms_fetch.py`, `firms_recent.py`: download and aggregate the NASA data.
- `geo_prep.py`, `harmonise.py`, `singapore.py`: build the harmonised record and the evidence.

## Team and tools

Solo entry by Yi Hang. Static site: HTML, CSS and JavaScript with no framework. Analysis: Python (pandas, NumPy, SciPy, Shapely).
AI assistance: Claude (Anthropic) was used for code, analysis design and drafting. Every figure on the site is computed from the data above by the scripts in this repository.
