"""How far to trust the wave model at each named break: its morning Hs against wave buoys.

Buoys come from NOAA's National Data Buoy Center (NDBC), which carries its own 3 m discus and
foam buoys, the Scripps CDIP Waveriders and the Great Lakes (GLOS) buoys, with no login. Each
is hand-matched in BUOYS to the named breaks within ~150 km on the same coast and exposure, and
in COASTS to the Lake Michigan and Lake Superior coasts of the regional pages (whose lake standard
is judged on total waves, which is what this compares everywhere).
Their standard meteorological records give WVHT (significant height, m), DPD (dominant period, s)
and MWD (direction, degrees), from four places on https://www.ndbc.noaa.gov/data/:

    historical/stdmet/{station}h{year}.txt.gz   complete years, quality-controlled
    stdmet/{Mon}/{station}{m}{year}.txt.gz      this year's finished months (m: 1-9, a, b, c)
    stdmet/{Mon}/{station}.txt                  the latest month, provisional
    realtime2/{station}.txt                     the last 45 days, provisional

Each morning that points.mornings() rates is paired with the buoy reading nearest the UTC instant
of the model step it used (MFWAM's 3-hourly step, or CDIP's hour), within MATCH. Over the matched
mornings, for the model's total Hs against the buoy's WVHT: bias (model minus buoy) and RMSE in
metres, Pearson r, and the share of mornings at 4 ft (1.22 m, also the lake standard's floor) or
more by each. On the mornings the buoy has 4 ft or more, bias4ft_m is the model's bias and hit4ft
the % the model also puts at 4 ft: the days that decide whether a morning rates at all, which a
small overall bias can hide. Each pair gets a trust label:

    good   |bias| < 0.2 m and r > 0.8
    fair   |bias| < 0.4 m and r > 0.6
    poor   anything else
    few    under MIN_N matched mornings (not judged)

A buoy is a point and the model a cell, often at a different depth and distance from shore, so
bias mixes model error with real differences in exposure; r is the cleaner test of whether the
model tracks the swell. In California, mornings() uses CDIP's MOP model at a 10-15 m deep site,
shoaled and island-shadowed, while most buoys sit in 20-500 m of water; MOP is also driven by
CDIP's own buoys, so those pairs are partly circular. California pairs therefore also carry
`offshore`: the same statistics for the raw Copernicus (MFWAM) cell that MOP replaced. Great Lakes
buoys are pulled out for the ice season, so their pairs cover roughly April to November.

    uv run python -m surfmap.buoys     # fetch() then validate(): writes data/buoys/validation.json
"""

import gzip
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from functools import lru_cache

import numpy as np
import requests
import xarray as xr

from . import cdip, config, points, rating

NDBC = "https://www.ndbc.noaa.gov/data"
OUT = config.DATA / "buoys"
RAW = OUT / "raw"
FIRST_YEAR = 2016                 # the morning series start on 31 December 2016
MATCH = np.timedelta64(60, "m")   # furthest a buoy reading may be from the model step
FT4 = 4 / rating.FT               # 1.22 m
MIN_N = 100                       # fewer matched mornings than this aren't judged
GOOD = (0.2, 0.8)                 # (largest |bias| m, smallest r) for "good"
FAIR = (0.4, 0.6)                 # and for "fair"
MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
HEADERS = {"User-Agent": "surfmap buoy validation (python-requests)"}

# (station, name, lat, lon, breaks): breaks are spots.SPOTS names on the same coast and exposure.
# Not used: 46206 La Perouse Bank (Tofino), 44258 Halifax Harbour (Lawrencetown), 62093 M4 Donegal
# Bay (Bundoran), 62107 Seven Stones (Fistral) and 22105 (Yangyang) are listed by NDBC but have no
# archived stdmet years; 46216 Goleta Point, 46223 Dana Point, 51204 Barbers Point, 45171 Granite
# Island and 41116/41118 off Sebastian Inlet ended in 2016-2018.
BUOYS = [
    # Hawaii: CDIP Waveriders off each break's own shore.
    ("51201", "Waimea Bay, Oahu (CDIP 106)", 21.671, -158.118, ("Pipeline", "Sunset Beach")),
    ("51208", "Hanalei, Kauai (CDIP 202)", 22.285, -159.574, ("Hanalei Bay",)),
    ("51205", "Pauwela, Maui (CDIP 187)", 21.018, -156.421, ("Honolua Bay",)),
    ("51211", "Pearl Harbor Entrance, Oahu (CDIP 233)", 21.297, -157.959, ("Waikiki",)),
    # California, north to south.
    ("46012", "Half Moon Bay, 24 nm SSW of San Francisco", 37.356, -122.881, ("Montara",)),
    ("46026", "San Francisco, 18 nm W", 37.750, -122.838, ("Montara", "Ocean Beach SF")),
    ("46237", "San Francisco Bar (CDIP 142)", 37.788, -122.634, ("Ocean Beach SF",)),
    ("46269", "Point Santa Cruz (CDIP 254)", 36.934, -122.034, ("Steamer Lane",)),
    ("46042", "Monterey, 27 nm WNW", 36.787, -122.408, ("Steamer Lane",)),
    ("46053", "East Santa Barbara, 12 nm SW", 34.246, -119.842, ("Rincon", "Ventura")),
    ("46217", "Anacapa Passage (CDIP 111)", 34.167, -119.435, ("Rincon", "Ventura")),
    ("46221", "Santa Monica Bay (CDIP 028)", 33.860, -118.641, ("Malibu",)),
    ("46268", "Topanga Nearshore (CDIP 103)", 34.022, -118.578, ("Malibu",)),
    ("46253", "San Pedro South (CDIP 213)", 33.576, -118.182, ("Huntington Beach",)),
    ("46222", "San Pedro (CDIP 092)", 33.614, -118.314, ("Huntington Beach",)),
    ("46277", "Green Beach Offshore (CDIP 271)", 33.336, -117.659, ("Trestles",)),
    ("46242", "Camp Pendleton Nearshore (CDIP 043)", 33.220, -117.439, ("Trestles",)),
    ("46224", "Oceanside Offshore (CDIP 045)", 33.178, -117.472, ("Trestles",)),
    ("46254", "Scripps Nearshore (CDIP 201)", 32.868, -117.267, ("Black's Beach",)),
    ("46225", "Torrey Pines Outer (CDIP 100)", 32.933, -117.391, ("Black's Beach",)),
    ("46232", "Point Loma South (CDIP 191)", 32.517, -117.425, ("Todos Santos",)),
    # Pacific Northwest.
    ("46278", "Tillamook Bay South Jetty (CDIP 270)", 45.559, -123.991, ("Short Sands",)),
    ("46029", "Columbia River Bar, 20 nm W", 46.148, -124.508, ("Short Sands",)),
    ("46211", "Grays Harbor (CDIP 036)", 46.857, -124.243, ("Westport",)),
    ("46041", "Cape Elizabeth, 45 nm NW of Aberdeen", 47.351, -124.741, ("Westport",)),
    ("46087", "Neah Bay, 6 nm N of Cape Flattery", 48.493, -124.727, ("Tofino",)),
    # Atlantic and Caribbean.
    ("41114", "Fort Pierce (CDIP 134)", 27.564, -80.215, ("Sebastian Inlet",)),
    ("41009", "Canaveral, 20 nm E", 28.508, -80.185, ("Sebastian Inlet",)),
    ("41120", "Cape Hatteras East (CDIP 250)", 35.258, -75.285, ("Cape Hatteras",)),
    ("41025", "Diamond Shoals", 35.026, -75.380, ("Cape Hatteras",)),
    ("44017", "Montauk Point, 23 nm SSW", 40.693, -72.049, ("Montauk",)),
    ("44097", "Block Island (CDIP 154)", 40.966, -71.123, ("Montauk",)),
    ("41115", "Rincon, Puerto Rico (CDIP 181)", 18.376, -67.280, ("Rincón",)),
    ("41121", "Arecibo, Puerto Rico (CDIP 249)", 18.491, -66.701, ("Rincón",)),
    # Western Pacific.
    ("52213", "North Point, Pohnpei (CDIP 178)", 7.081, 158.244, ("P-Pass",)),
    # Great Lakes.
    ("45218", "Sheboygan, Wisconsin Shipwreck Coast", 43.731, -87.624, ("Sheboygan",)),
    ("45210", "Rawley Point East (CDIP 269)", 44.055, -87.050, ("Sheboygan",)),
    ("45013", "Atwater Park, Milwaukee", 43.098, -87.850, ("Sheboygan",)),
    ("45161", "Muskegon", 43.184, -86.358, ("Grand Haven",)),
    ("45029", "Holland", 42.900, -86.272, ("Grand Haven",)),
    ("45027", "McQuade Harbor Nearshore", 46.860, -91.930, ("Stoney Point",)),
    ("45028", "Western Lake Superior", 46.814, -91.829, ("Stoney Point",)),
    ("45173", "Munising", 46.573, -86.572, ("Marquette",)),
    ("45211", "Grand Island North (CDIP 268)", 46.973, -86.568, ("Marquette",)),
    ("45004", "East Superior, 70 nm NE of Marquette", 47.583, -86.586, ("Marquette",)),
    ("45176", "Cleveland Intake Crib", 41.550, -81.765, ("Edgewater",)),
    ("45205", "Edgewater Beach", 41.501, -81.748, ("Edgewater",)),
    # Great Lakes buoys used only for the regional coasts (COASTS).
    ("45198", "Chicago", 41.892, -87.563, ()),
    ("45177", "Ohio Street Beach, Chicago", 41.894, -87.613, ()),
    ("45174", "Wilmette", 42.135, -87.655, ()),
    ("45187", "Winthrop Harbor", 42.491, -87.779, ()),
    ("45199", "Salmon Unlimited, Racine", 42.702, -87.647, ()),
    ("45007", "South Michigan, 43 nm ESE of Milwaukee", 42.674, -87.026, ()),
    ("45170", "Michigan City", 41.755, -86.968, ()),
    ("45026", "Cook Nuclear Plant, Stevensville", 41.982, -86.619, ()),
    ("45168", "South Haven", 42.397, -86.331, ()),
    ("45024", "Ludington", 43.980, -86.560, ()),
    ("45183", "Sleeping Bear Dunes", 44.982, -85.831, ()),
    ("45002", "North Michigan, mid-lake", 45.344, -86.411, ()),
    ("45006", "West Superior, 30 nm NE of Outer Island", 47.335, -89.793, ()),
    ("45216", "Ontonagon", 46.907, -89.354, ()),
    ("45023", "North Entry, Keweenaw Peninsula", 47.270, -88.607, ()),
    ("45172", "Grand Marais, Michigan", 46.740, -85.980, ()),
    ("45001", "Mid Superior, 60 nm NNE of Hancock", 48.061, -87.793, ()),
]

# {regional page: {coast: buoys}} for regions.REGIONS coasts, nearest first, then a mid-lake buoy
# where it lies in the coast's fetch. Not used: 45136 Slate Island (Terrace Bay) and 45217
# Wisconsin Point (Park Point) have no archived years, 45219 Two Harbors no wave heights, and
# 45025 South Entry faces into Keweenaw Bay.
COASTS = {
    "Lake Michigan": {
        "Whiting": ("45198", "45170", "45007"), "Gary": ("45170", "45198", "45007"), "Portage": ("45170", "45007"),
        "Chicago": ("45177", "45198", "45174"), "Racine": ("45199", "45187", "45007"), "Milwaukee": ("45013", "45007"),
        "Sheboygan": ("45218", "45210", "45013"), "New Buffalo": ("45170", "45026"), "St. Joseph": ("45026", "45168"),
        "Holland": ("45029", "45168", "45007"), "Grand Haven": ("45161", "45029", "45007"),
        "Ludington": ("45024", "45210"), "Frankfort": ("45183", "45024", "45002"),
    },
    "Lake Superior": {
        "Park Point": ("45027", "45028"), "Stoney Point": ("45028", "45027"), "Grand Marais MN": ("45006",),
        "Little Girl's Point": ("45216", "45006"), "Ontonagon": ("45216", "45006", "45023"),
        "Marquette": ("45173", "45211", "45004"), "Au Train": ("45173", "45211", "45004"),
        "Grand Marais MI": ("45172", "45211", "45004"), "Whitefish Point": ("45172", "45004"),
        "Agawa Bay": ("45004",), "Terrace Bay": ("45001",),
    },
}


def _get(url):
    """The body at url, or None where NDBC has no such file; other failures retry with backoff."""
    for attempt in range(5):
        try:
            r = requests.get(url, headers=HEADERS, timeout=60)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.content
        except requests.RequestException:
            if attempt == 4:
                raise
            time.sleep(5 * 2**attempt)


def _save(path, url, max_age_days=None):
    """Cache url at path (empty where NDBC has none); re-download once older than max_age_days."""
    if path.exists() and (max_age_days is None or time.time() - path.stat().st_mtime < max_age_days * 86400):
        return
    body = _get(url)
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_bytes(body or b"")
    tmp.rename(path)


@lru_cache(maxsize=None)
def _month_files(mon):
    """File names in NDBC's stdmet/{mon}/ directory: this year's months, and the latest provisional one."""
    html = (_get(f"{NDBC}/stdmet/{mon}/") or b"").decode("ascii", "replace")
    return frozenset(re.findall(r'href="([0-9a-z]+(?:[1-9abc]\d{4}\.txt\.gz|\.txt))"', html))


def fetch_station(station):
    """Every year of one buoy's record from FIRST_YEAR on: whole years, then this year's months."""
    folder = RAW / station
    folder.mkdir(parents=True, exist_ok=True)
    this_year = date.today().year
    have = set()
    for year in range(FIRST_YEAR, this_year):
        path = folder / f"{station}h{year}.txt.gz"
        recent = year >= this_year - 1 and path.exists() and not path.stat().st_size
        _save(path, f"{NDBC}/historical/stdmet/{path.name}", 7 if recent else None)  # last year's may come late
        if path.stat().st_size:
            have.add(year)
    for k, mon in enumerate(MONTHS):
        files = _month_files(mon)
        name = f"{station}{'123456789abc'[k]}"
        for year in range(FIRST_YEAR, this_year + 1):
            if f"{name}{year}.txt.gz" in files and year not in have:
                _save(folder / f"{name}{year}.txt.gz", f"{NDBC}/stdmet/{mon}/{name}{year}.txt.gz")
        if f"{station}.txt" in files:
            _save(folder / f"{station}_{mon}.txt", f"{NDBC}/stdmet/{mon}/{station}.txt", 1)
    _save(folder / f"{station}_realtime.txt", f"{NDBC}/realtime2/{station}.txt", 1)


def fetch(workers=4):
    """Download every buoy in BUOYS into data/buoys/raw/, a few stations at a time."""
    for mon in MONTHS:  # list the month directories once, before the threads share them
        _month_files(mon)
    stations = [b[0] for b in BUOYS]
    with ThreadPoolExecutor(workers) as pool:
        for station, err in zip(stations, pool.map(_try_fetch, stations)):
            print(f"  {station}: {err or 'ok'}", flush=True)


def _try_fetch(station):
    try:
        fetch_station(station)
    except Exception as e:  # keep going; a re-run picks up what failed
        return f"failed: {e}"


def _read(path):
    """(time, wvht, dpd, mwd) from one NDBC stdmet file, with missing values as NaN."""
    raw = path.read_bytes()
    if not raw:
        return None
    lines = (gzip.decompress(raw) if path.suffix == ".gz" else raw).decode("ascii", "replace").splitlines()
    head = lines[0].lstrip("#").split()
    body = [line for line in lines if line.strip() and not line.startswith("#")]
    if any(len(line.split()) != len(head) for line in body):  # a truncated line: drop it
        body = [line for line in body if len(line.split()) == len(head)]
    if not body:
        return None
    a = np.array(" ".join(body).replace("MM", "nan").split(), dtype=float).reshape(len(body), len(head))
    col = {h: a[:, i] for i, h in enumerate(head)}
    year = col["YY"] if "YY" in col else col["YYYY"]
    month = ((year - 1970) * 12 + col["MM"] - 1).astype(int).astype("datetime64[M]")
    minutes = ((col["DD"] - 1) * 1440 + col["hh"] * 60 + col.get("mm", 0 * year)).astype(int)
    t = month.astype("datetime64[m]") + minutes.astype("timedelta64[m]")
    hs, dpd, mwd = col["WVHT"], col["DPD"], col["MWD"]
    return t, np.where((hs >= 0) & (hs < 25), hs, np.nan), np.where(dpd < 99, dpd, np.nan), np.where(mwd <= 360, mwd, np.nan)


@lru_cache(maxsize=None)
def series(station):
    """One buoy's wave record in time order: dict(t, hs, dpd, mwd), readings without WVHT dropped.
    Where files overlap (a provisional month and its final version) the final one wins."""
    folder = RAW / station
    final = sorted(folder.glob("*.txt.gz"))
    provisional = sorted(folder.glob("*_???.txt")) + sorted(folder.glob("*_realtime.txt"))
    parts = [p for p in map(_read, final + provisional) if p is not None]
    if not parts:
        return None
    t, hs, dpd, mwd = (np.concatenate(x) for x in zip(*parts))
    ok = np.isfinite(hs)
    t, hs, dpd, mwd = t[ok], hs[ok], dpd[ok], mwd[ok]
    order = np.argsort(t, kind="stable")
    keep = order[np.unique(t[order], return_index=True)[1]]
    return dict(t=t[keep], hs=hs[keep], dpd=dpd[keep], mwd=mwd[keep])


def _samples(lat, lon):
    """(day, utc time, hs) of the Copernicus step behind each of points.mornings()' mornings,
    and the utc time of each morning as mornings() reports it (CDIP's hour wherever it took over)."""
    with points._open([points._waves_path(lat, lon, True), points._waves_path(lat, lon)]) as w:
        idx, days = points._local_mornings(w.time.values, lon, 3)
        t = w.time.values[idx].astype("datetime64[m]")
        hs = w["VHM0"].values[idx].astype(float)
    used = t.copy()
    site = cdip.site_for(lat, lon)
    if site and cdip._path(site).exists():  # the same swap as points.mornings()
        with xr.open_dataset(cdip._path(site)) as d:
            ct, chs = d.waveTime.values.astype("datetime64[m]"), d.waveHs.values
        cday = (ct + np.timedelta64(int(round(lon / 15 * 60)), "m")).astype("datetime64[D]")
        pos = np.minimum(np.searchsorted(cday, days), len(cday) - 1)
        hit = (cday[pos] == days) & np.isfinite(chs[pos])
        used = np.where(hit, ct[pos], t)
    return days, t, hs, used


def _nearest(t, tb):
    """Index into the sorted times tb of the reading nearest each of t, or -1 if none within MATCH."""
    k = np.clip(np.searchsorted(tb, t), 1, len(tb) - 1)
    k -= (t - tb[k - 1]) < (tb[k] - t)
    return np.where(np.abs(tb[k] - t) <= MATCH, k, -1)


def trust(n, bias, r):
    """The label in the module docstring, from matched mornings, bias (m) and correlation."""
    if n < MIN_N:
        return "few"
    if abs(bias) < GOOD[0] and r > GOOD[1]:
        return "good"
    if abs(bias) < FAIR[0] and r > FAIR[1]:
        return "fair"
    return "poor"


def compare(model_hs, t, buoy):
    """Statistics of model Hs at times t against the buoy record (series()): matched mornings,
    bias and RMSE (m), r, % of mornings at 4 ft or more by each, the trust label, first/last day."""
    k = _nearest(t, buoy["t"])
    ok = (k >= 0) & np.isfinite(model_hs)
    m, b = model_hs[ok], buoy["hs"][k[ok]]
    n = int(ok.sum())
    if n < 3:
        return dict(n=n, bias_m=None, rmse_m=None, r=None, pct4ft_buoy=None, pct4ft_model=None,
                    hs_buoy_m=None, bias4ft_m=None, hit4ft=None, trust="few", first=None, last=None)
    d, big = m - b, b >= FT4
    out = dict(bias_m=float(d.mean()), rmse_m=float(np.sqrt(np.mean(d**2))), r=float(np.corrcoef(m, b)[0, 1]),
               pct4ft_buoy=float(100 * np.mean(big)), pct4ft_model=float(100 * np.mean(m >= FT4)),
               hs_buoy_m=float(b.mean()), bias4ft_m=float(d[big].mean()) if big.any() else None,
               hit4ft=float(100 * np.mean(m[big] >= FT4)) if big.any() else None)
    out = {key: v if v is None else round(v, 3) for key, v in out.items()}
    out["trust"] = trust(n, out["bias_m"], out["r"])
    days = t[ok].astype("datetime64[D]")
    return dict(n=n, **out, first=str(days.min()), last=str(days.max()))


def _pairs():
    """[(kind, name, region, lat, lon, buoy)] for the breaks in BUOYS (in spots.SPOTS order), then
    the coasts in COASTS, each point's buoys nearest first."""
    pts = {(kind, name, region): (lat, lon) for kind, name, region, lat, lon, _ in points.all_points()}
    out = []
    for (kind, name, region), (lat, lon) in pts.items():
        if kind == "break":
            near = [b for b in BUOYS if name in b[4]]
        elif kind == "coast" and name in COASTS.get(region, {}):
            near = [b for b in BUOYS if b[0] in COASTS[region][name]]
        else:
            continue
        near.sort(key=lambda b: cdip._km(lat, lon, b[2], b[3]))
        out += [(kind, name, region, lat, lon, b) for b in near]
    return out


def validate():
    """Compare every break-buoy and coast-buoy pair and write data/buoys/validation.json; returns the list."""
    out = []
    for kind, name, region, lat, lon, (station, buoy_name, blat, blon, _) in _pairs():
        m, buoy = points.mornings(lat, lon), series(station)
        if m is None or buoy is None or len(buoy["t"]) < 2:
            print(f"  {name} / {station}: no {'model' if m is None else 'buoy'} data", flush=True)
            continue
        days, t_off, hs_off, t_used = _samples(lat, lon)
        assert np.array_equal(days, m["day"]), name
        row = {"break": name, "kind": kind, "region": region, "buoy": station, "buoy_name": buoy_name, "buoy_lat": blat,
               "buoy_lon": blon, "km": round(float(cdip._km(lat, lon, blat, blon)), 1), "source": m["source"]}
        row.update(compare(m["hs"], t_used, buoy))
        if m["source"].startswith("CDIP"):
            row["offshore"] = dict(source="MFWAM ~8 km cell (before the CDIP swap)", **compare(hs_off, t_off, buoy))
        out.append(row)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "validation.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    _print(out)
    return out


def _print(rows):
    f = lambda v, fmt: "   -" if v is None else format(v, fmt)
    print(f"{'break':18s} {'buoy':6s} {'km':>5s} {'n':>5s} {'bias':>6s} {'rmse':>5s} {'r':>5s} {'4ft b/m %':>10s} trust"
          "   | offshore bias, r, trust")
    for k, r in enumerate(rows):
        if r["kind"] == "coast" and (k == 0 or rows[k - 1]["region"] != r["region"]):
            print(f"-- {r['region']} coasts")
        line = (f"{r['break']:18s} {r['buoy']:6s} {r['km']:5.0f} {r['n']:5d} {f(r['bias_m'], '+6.2f')} "
                f"{f(r['rmse_m'], '5.2f')} {f(r['r'], '5.2f')} {f(r['pct4ft_buoy'], '4.0f')}/{f(r['pct4ft_model'], '<4.0f')} "
                f"{r['trust']:5s}")
        if "offshore" in r:
            o = r["offshore"]
            line += f" | {f(o['bias_m'], '+6.2f')} {f(o['r'], '5.2f')} {o['trust']}"
        print(line)


if __name__ == "__main__":
    fetch()
    validate()
