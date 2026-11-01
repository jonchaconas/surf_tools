"""Open-Meteo marine (Meteo-France MFWAM) for the basins WW3 global leaves blank.

Open-Meteo's free tier is metered in weighted calls (each location counts once
per 14 days requested; limits are 600/min, 5,000/h and 10,000/day), so these
basins are sampled on a coarse grid of water points rather than every model cell.
"""

import json
import time
from functools import lru_cache

import cartopy.io.shapereader as shp
import numpy as np
import requests
from shapely import contains_xy
from shapely.geometry import box
from shapely.ops import unary_union

from . import config, ww3

API = "https://marine-api.open-meteo.com/v1/marine"
# Up to 10 variables cost the same quota, so take the total-sea fields too.
HOURLY = ("swell_wave_height", "swell_wave_period", "swell_wave_peak_period", "swell_wave_direction",
          "wind_wave_height", "wave_height", "wave_period", "wave_direction")
OUT = config.DATA / "inland"
CHUNK = 20         # locations per request
PER_MINUTE = 450   # weighted calls, under the free tier's 600
PER_HOUR = 4500    # under the free tier's 5,000; the 10,000/day cap is handled by waiting on 429s

GREAT_LAKES = ("Lake Superior", "Lake Michigan", "Lake Huron", "Lake Erie", "Lake Ontario")

# name: (grid spacing deg, lat0, lat1, lon0, lon1). Earlier boxes take precedence where they overlap.
BASINS = {
    "great_lakes": (0.75, 41.0, 49.5, -92.5, -76.0),
    "black_sea": (1.5, 40.9, 47.5, 27.0, 42.0),
    "caspian": (1.5, 36.5, 47.5, 46.5, 55.5),
    "mediterranean": (2.0, 30.0, 46.0, -5.6, 36.5),
}
# North of 42N the Mediterranean box would take in the Bay of Biscay, so it starts at 0E there.
ATLANTIC_CUT = box(-5.6, 42.0, 0.0, 46.0)


@lru_cache
def geometries():
    """Water polygon for each basin, from Natural Earth 1:50m."""
    lakes = shp.Reader(shp.natural_earth("50m", "physical", "lakes")).records()
    lakes = unary_union([r.geometry for r in lakes if r.attributes["name"] in GREAT_LAKES])
    ocean = unary_union(list(shp.Reader(shp.natural_earth("50m", "physical", "ocean")).geometries()))
    out, taken = {}, None
    for name, (_, la0, la1, lo0, lo1) in BASINS.items():
        area = box(lo0, la0, lo1, la1)
        geom = (lakes if name == "great_lakes" else ocean).intersection(area)
        if name == "mediterranean":
            geom = geom.difference(ATLANTIC_CUT)
        out[name] = geom if taken is None else geom.difference(taken)
        taken = area if taken is None else taken.union(area)
    return out


def points(name):
    """Grid points (lat, lon) inside a basin's water polygon."""
    d, la0, la1, lo0, lo1 = BASINS[name]
    la, lo = np.meshgrid(np.arange(la0 + d / 2, la1, d), np.arange(lo0 + d / 2, lo1, d), indexing="ij")
    keep = contains_xy(geometries()[name], lo, la)
    return np.round(la[keep], 3), np.round(lo[keep], 3)


def _window(year, month):
    """Date range matching the days WW3 has, so both sources cover the same mornings."""
    stamps = ww3.available_hours(year, month)
    return stamps[0].date(), stamps[-1].date()


class _Throttle:
    """Keeps this process's weighted calls under per-minute and per-hour budgets."""

    def __init__(self, per_minute, per_hour):
        self.limits, self.log = ((60, per_minute), (3600, per_hour)), []

    def wait(self, weight):
        while True:
            now = time.monotonic()
            self.log = [(t, w) for t, w in self.log if now - t < 3600]
            if all(sum(w for t, w in self.log if now - t < span) + weight <= cap
                   for span, cap in self.limits):
                break
            time.sleep(5)
        self.log.append((time.monotonic(), weight))


def _request(lats, lons, d0, d1, tries=200):
    params = dict(
        latitude=",".join(map(str, lats)),
        longitude=",".join(map(str, lons)),
        start_date=d0.isoformat(),
        end_date=d1.isoformat(),
        hourly=",".join(HOURLY),
        timezone="GMT",
    )
    for _ in range(tries):
        try:
            r = requests.get(API, params=params, timeout=180)
        except (requests.Timeout, requests.ConnectionError) as e:
            print(f"  Open-Meteo {type(e).__name__}; retrying", flush=True)
            time.sleep(30)
            continue
        if r.status_code == 429 or not r.content:
            reason = r.text[:160]
            wait = 3600 if "daily" in reason.lower() else 600 if "hourly" in reason.lower() else 65
            print(f"  Open-Meteo throttled ({r.status_code}: {reason}); waiting {wait}s", flush=True)
            time.sleep(wait)
            continue
        r.raise_for_status()
        data = r.json()
        return data if isinstance(data, list) else [data]
    raise RuntimeError("Open-Meteo kept refusing; rerun fetch later, finished chunks are cached")


def _files(name, year, month):
    n = len(points(name)[0])
    return [OUT / f"{name}_{year}-{month:02d}_{k:02d}.json" for k in range(-(-n // CHUNK))]


def complete(name, month):
    """True once every season of this month is cached for the basin."""
    return all(p.exists() for y in config.seasons(month) for p in _files(name, y, month))


def fetch():
    """Fills whole months in page order, so a partial run still yields finished pages."""
    OUT.mkdir(parents=True, exist_ok=True)
    throttle, total = _Throttle(PER_MINUTE, PER_HOUR), 0.0
    pts = {name: points(name) for name in BASINS}
    for month in config.MONTHS:
        for year in config.seasons(month):
            d0, d1 = _window(year, month)
            per_loc = max(1.0, ((d1 - d0).days + 1) / 14)
            for name, (lats, lons) in pts.items():
                for k, path in enumerate(_files(name, year, month)):
                    if path.exists():
                        continue
                    la, lo = lats[k * CHUNK:(k + 1) * CHUNK], lons[k * CHUNK:(k + 1) * CHUNK]
                    throttle.wait(per_loc * len(la))
                    rows = _request(la, lo, d0, d1)
                    for row, a, b in zip(rows, la, lo):
                        row["requested"] = [float(a), float(b)]
                    path.write_text(json.dumps(rows))
                    total += per_loc * len(la)
                    print(f"  {path.name} ({total:.0f} weighted calls so far)", flush=True)


def morning(row, keep):
    """One location's hourly response cut to its ~config.LOCAL_HOUR sample on the days in keep."""
    lon = row["requested"][1]
    t = np.array(row["hourly"]["time"], dtype="datetime64[h]")
    day = t.astype("datetime64[D]")
    h = round(config.LOCAL_HOUR - lon / 15) % 24
    sel = ((t - day).astype(int) == h) & np.isin(day, keep)
    return {v: np.array(row["hourly"][v], dtype=float)[sel] for v in HOURLY}


def load(month):
    """Morning samples per basin: {name: dict(lat, lon, hs, tp, dir, ws)}, arrays shaped (point, day).

    Open-Meteo has no swell peak period over the Great Lakes, so there it is estimated
    from the mean swell period using the median peak/mean ratio of the other basins;
    that ratio is returned as tp_ratio on the estimated basins.
    """
    keep = np.array(ww3.days(month), dtype="datetime64[D]")
    raw = {}
    for name in BASINS:
        if not complete(name, month):
            continue
        series = {}
        for path in sorted(OUT.glob(f"{name}_*-{month:02d}_*.json")):
            for row in json.loads(path.read_text()):
                series.setdefault(tuple(row["requested"]), []).append(morning(row, keep))
        if series:
            pts = list(series)
            raw[name] = dict(lat=np.array([p[0] for p in pts]), lon=np.array([p[1] for p in pts]),
                             **{v: np.stack([np.concatenate([s[v] for s in series[p]]) for p in pts])
                                for v in HOURLY})

    if not raw:
        return {}
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.nanmedian(np.concatenate(
            [(b["swell_wave_peak_period"] / b["swell_wave_period"]).ravel() for b in raw.values()]))
    out = {}
    for name, b in raw.items():
        tp, est = b["swell_wave_peak_period"], np.isnan(b["swell_wave_peak_period"]).all()
        out[name] = dict(lat=b["lat"], lon=b["lon"], hs=b["swell_wave_height"],
                         tp=b["swell_wave_period"] * ratio if est else tp,
                         dir=b["swell_wave_direction"], ws=b["wind_wave_height"],
                         tp_ratio=ratio if est else None)
    return out
