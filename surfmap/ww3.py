"""WAVEWATCH III global (PacIOOS ERDDAP), one morning snapshot per day at 1 degree.

The model has no data for the Mediterranean, Black Sea, Caspian or Great Lakes;
inland.py covers those.
"""

import calendar
import os
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

import numpy as np
import requests
import xarray as xr

from . import config

ERDDAP = "https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global"
VARS = ("shgt", "sper", "sdir", "whgt")
BAND = 30    # degrees of longitude per request; each band is sampled at its own UTC hour
STRIDE = 2   # native 0.5 degree grid -> 1 degree
OUT = config.DATA / "ww3"


def band_hour(lon0):
    """UTC hour nearest config.LOCAL_HOUR solar time at the band's centre longitude."""
    centre = lon0 + BAND / 2
    if centre > 180:
        centre -= 360
    return round(config.LOCAL_HOUR - centre / 15) % 24


def _quote(q):
    return q.replace("[", "%5B").replace("]", "%5D")


def _download(url, path, tries=4):
    for attempt in range(tries):
        try:
            r = requests.get(url, timeout=900)
            r.raise_for_status()
            tmp = path.with_name(f"{path.name}.{os.getpid()}.part")
            tmp.write_bytes(r.content)
            tmp.rename(path)
            return
        except requests.RequestException:
            if attempt == tries - 1:
                raise
            time.sleep(15 * (attempt + 1))


def available_hours(year, month):
    """Timestamps the model has for one month; the archive has a few outages."""
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"time_{year}-{month:02d}.csv"
    if not path.exists():
        last = calendar.monthrange(year, month)[1]
        q = f"time[({year}-{month:02d}-01T00:00:00Z):1:({year}-{month:02d}-{last}T23:00:00Z)]"
        _download(f"{ERDDAP}.csv?{_quote(q)}", path)
    stamps = (datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ") for s in path.read_text().split()[2:])
    # ERDDAP snaps the range to the nearest stamps, which can fall outside the month.
    return sorted(t for t in stamps if (t.year, t.month) == (year, month))


def _segments(stamps):
    """Split hourly stamps into gap-free runs, so a 24-step index stride stays on one hour."""
    runs = [[stamps[0]]]
    for t in stamps[1:]:
        if t - runs[-1][-1] == timedelta(hours=1):
            runs[-1].append(t)
        else:
            runs.append([t])
    return runs


def _jobs():
    for year in config.YEARS:
        for month in config.MONTHS:
            for seg in _segments(available_hours(year, month)):
                for lon0 in range(0, 360, BAND):
                    h = band_hour(lon0)
                    picks = [t for t in seg if t.hour == h]
                    if picks:
                        yield year, month, lon0, h, picks[0], picks[-1]


def _path(year, month, lon0, t0):
    return OUT / f"{year}-{month:02d}_lon{lon0:03d}_d{t0.day:02d}.nc"


def _url(lon0, t0, t1):
    sel = (
        f"[({t0:%Y-%m-%dT%H}:00:00Z):24:({t1:%Y-%m-%dT%H}:00:00Z)][(0.0)]"
        f"[(-77.5):{STRIDE}:(77.5)][({lon0}):{STRIDE}:({lon0 + BAND - 0.5})]"
    )
    return f"{ERDDAP}.nc?" + ",".join(_quote(v + sel) for v in VARS)


def fetch(workers=3):
    jobs = [j for j in _jobs() if not _path(j[0], j[1], j[2], j[4]).exists()]
    print(f"WW3: {len(jobs)} files to download", flush=True)

    def one(job):
        year, month, lon0, h, t0, t1 = job
        path = _path(year, month, lon0, t0)
        _download(_url(lon0, t0, t1), path)
        with xr.open_dataset(path) as ds:
            hours = set(ds.time.dt.hour.values.tolist())
        if hours != {h}:
            path.unlink()
            raise RuntimeError(f"{path.name}: expected hour {h}, got {sorted(hours)}")
        print(f"  {path.name}", flush=True)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, jobs))


def load(month):
    """Daily morning samples for one calendar month, all years: dims (day, lat, lon).

    Longitudes are converted to -180..180. Bands can differ in which days exist
    (outages fall on different UTC hours), so missing band-days are NaN.
    """
    bands = []
    for lon0 in range(0, 360, BAND):
        parts = []
        for year in config.YEARS:
            for path in sorted(OUT.glob(f"{year}-{month:02d}_lon{lon0:03d}_d*.nc")):
                ds = xr.open_dataset(path).squeeze("depth", drop=True)
                ds = ds.assign_coords(day=ds.time.dt.floor("D")).swap_dims(time="day").drop_vars("time")
                parts.append(ds[list(VARS)].load())
        if parts:
            bands.append(xr.concat(parts, "day"))
    ds = xr.concat(bands, "longitude", join="outer")
    ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180).sortby("longitude")
    return ds.rename(latitude="lat", longitude="lon")


def days(month):
    """Every UTC date in the month that has at least one model hour."""
    return sorted({t.date() for y in config.YEARS for t in available_hours(y, month)})


def spot_cell(ds, lat, lon, max_km=250):
    """Nearest ocean cell to a coastal point, or None when nothing is within max_km."""
    valid = np.isfinite(ds.shgt).mean("day") > 0.5
    la, lo = np.meshgrid(ds.lat.values, ds.lon.values, indexing="ij")
    dlon = (lo - lon + 180) % 360 - 180
    km = 111.2 * np.hypot(la - lat, dlon * np.cos(np.radians(lat)))
    km = np.where(valid.values, km, np.inf)
    i, j = np.unravel_index(np.argmin(km), km.shape)
    return None if km[i, j] > max_km else (i, j)
