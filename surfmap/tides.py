"""Tidal range at each break, from EOT20 (DGFI-TUM's global ocean tide model, open access).

EOT20 gives each tidal constituent's amplitude and phase on a 1/8 degree grid. At each point
the five main constituents are summed into a year of hourly tide, and the daily ranges are
summarised as a typical day (median) and a big-tide day (90th percentile). Astronomical
start phases are ignored: they shift when spring tides fall but not how big days get over a
year. Download: https://www.seanoe.org/data/00683/79489/ (ocean_tides/*_ocean_eot20.nc).
"""

import json
from functools import lru_cache

import numpy as np
import xarray as xr

from . import config

DIR = config.DATA / "eot20" / "ocean_tides"
CACHE = config.DATA / "eot20" / "ranges.json"
SPEEDS = {"M2": 28.9841042, "S2": 30.0, "N2": 28.4397295, "K1": 15.0410686, "O1": 13.9430356}  # degrees per hour
SEARCH = 0.5  # degrees: how far to look for the nearest cell with tide data


@lru_cache
def _grids():
    return {c: xr.open_dataset(DIR / f"{c}_ocean_eot20.nc") for c in SPEEDS}


def _constituents(lat, lon):
    """{constituent: (amplitude m, phase deg)} at the nearest cell with data, or None."""
    lon = lon % 360
    out = {}
    for c, ds in _grids().items():
        box = ds.sel(lat=slice(lat - SEARCH, lat + SEARCH), lon=slice(lon - SEARCH, lon + SEARCH))
        amp = box.amplitude.values
        if not np.isfinite(amp).any():
            return None
        la, lo = np.meshgrid(box.lat.values, box.lon.values, indexing="ij")
        km = np.where(np.isfinite(amp), np.hypot(la - lat, (lo - lon) * np.cos(np.radians(lat))), np.inf)
        i, j = np.unravel_index(np.argmin(km), km.shape)
        out[c] = (float(amp[i, j]) / 100, float(box.phase.values[i, j]))
    return out


def _ranges(cons):
    """(typical, big-tide) daily range in metres over a year of hourly tide."""
    t = np.arange(24 * 365)
    h = sum(a * np.cos(np.radians(SPEEDS[c] * t - ph)) for c, (a, ph) in cons.items())
    daily = h.reshape(365, 24)
    r = daily.max(1) - daily.min(1)
    return float(np.median(r)), float(np.percentile(r, 90))


def _points():
    from .regions import REGIONS
    from .spots import SPOTS
    from .stretches import STRETCHES

    coasts = [c for spec in REGIONS.values() for c in spec["coasts"]]
    return {name: (lat, lon) for name, _, lat, lon, *_ in SPOTS + STRETCHES + coasts}


@lru_cache
def ranges():
    """{name: (typical m, big-tide m)} for every break, stretch and regional coast (cached on disk)."""
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    points = _points()
    missing = [n for n in points if n not in cache]
    if missing and DIR.exists():
        for name in missing:
            cons = _constituents(*points[name])
            cache[name] = _ranges(cons) if cons else None
        CACHE.write_text(json.dumps(cache, indent=0))
    return {n: tuple(v) if v else None for n, v in cache.items()}


def text(name, ft):
    """'4–7 ft' (typical to big-tide range), or '' if unknown."""
    r = ranges().get(name)
    return f"{r[0] * ft:.0f}–{r[1] * ft:.0f} ft" if r else ""
