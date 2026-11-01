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


# Equilibrium arguments (Schureman): V0 = a*T + b*s + c*h + d*p + e degrees, with T = 15 * UT hours,
# s/h/p the mean longitudes of moon, sun and lunar perigee.
_ARGS = {"M2": (2, -2, 2, 0, 0), "S2": (2, 0, 0, 0, 0), "N2": (2, -3, 2, 1, 0), "K1": (1, 0, 1, 0, 90), "O1": (1, -2, 1, 0, -90)}


def _nodal(c, n):
    """Nodal amplitude factor f and phase u (degrees) for the 18.6-year lunar cycle, N in degrees."""
    cn, sn = np.cos(np.radians(n)), np.sin(np.radians(n))
    if c in ("M2", "N2"):
        return 1 - 0.037 * cn, -2.1 * sn
    if c == "K1":
        return 1.006 + 0.115 * cn, -8.9 * sn
    if c == "O1":
        return 1.009 + 0.187 * cn, 10.8 * sn
    return 1.0, 0.0


def height(lat, lon, times):
    """Tide height (m, about mean sea level) at a point at the given UTC times (datetime64), from the
    five EOT20 constituents with astronomical arguments and nodal corrections; None without data."""
    cons = _constituents(lat, lon)
    if cons is None:
        return None
    t = np.asarray(times, dtype="datetime64[s]")
    d = (t - np.datetime64("2000-01-01T12:00:00")) / np.timedelta64(1, "D")  # days from J2000
    hours = (t - t.astype("datetime64[D]")) / np.timedelta64(1, "h")
    big_t = 15 * hours  # checked against NOAA's predictions for San Francisco (r = 0.98)
    s, h, p, n = 218.3165 + 13.17639648 * d, 280.4661 + 0.98564736 * d, 83.3535 + 0.11140353 * d, 125.0445 - 0.05295377 * d
    out = 0
    for c, (amp, phase) in cons.items():
        a, b, cc, dd, e = _ARGS[c]
        f, u = _nodal(c, n)
        out = out + f * amp * np.cos(np.radians(a * big_t + b * s + cc * h + dd * p + e + u - phase))
    return out


def text(name, ft):
    """'4–7 ft' (typical to big-tide range), or '' if unknown."""
    r = ranges().get(name)
    return f"{r[0] * ft:.0f}–{r[1] * ft:.0f} ft" if r else ""
