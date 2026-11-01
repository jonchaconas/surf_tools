"""WAVEWATCH III global (PacIOOS ERDDAP), one morning snapshot per day at 1 degree.

The model has no data for the Mediterranean, Black Sea, Caspian or Great Lakes;
inland.py covers those. In DETAIL boxes it is also fetched at its native 0.5 degrees,
because 1 degree sampling misses most of the shallow seas there.
"""

import calendar
import os
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

import requests
import xarray as xr

from . import config

ERDDAP = "https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global"
ARCHIVE_START = datetime(2017, 1, 1, 21)  # first stamp in the archive (config.FIRST_YEAR)
VARS = ("shgt", "sper", "sdir", "whgt")
BAND = 30    # degrees of longitude per request; each band is sampled at its own UTC hour
STRIDE = 2   # native 0.5 degree grid -> 1 degree
OUT = config.DATA / "ww3"

# name: (lat0, lat1, lon0, lon1), longitudes 0-360 east.
DETAIL = {
    "east_asia": (-10.0, 42.0, 98.0, 131.0),  # Gulfs of Thailand and Tonkin, China's seas, Java Sea
}


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
        first = max(datetime(year, month, 1), ARCHIVE_START)  # ERDDAP 404s on a start before the archive's
        q = f"time[({first:%Y-%m-%dT%H:%M:%SZ}):1:({year}-{month:02d}-{last}T23:00:00Z)]"
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


def _path(year, month, lon0, t0, region=None):
    tag = f"_{region}" if region else ""
    return OUT / f"{year}-{month:02d}{tag}_lon{lon0:03d}_d{t0.day:02d}.nc"


def _url(t0, t1, lat0, lat1, lon0, lon1, stride):
    sel = (
        f"[({t0:%Y-%m-%dT%H}:00:00Z):24:({t1:%Y-%m-%dT%H}:00:00Z)][(0.0)]"
        f"[({lat0}):{stride}:({lat1})][({lon0}):{stride}:({lon1})]"
    )
    return f"{ERDDAP}.nc?" + ",".join(_quote(v + sel) for v in VARS)


def _jobs():
    """(path, url, UTC hour) for every file the configured months need."""
    for month in config.MONTHS:
        for year in config.seasons(month):
            for seg in _segments(available_hours(year, month)):
                for lon0 in range(0, 360, BAND):
                    h = band_hour(lon0)
                    picks = [t for t in seg if t.hour == h]
                    if not picks:
                        continue
                    t0, t1 = picks[0], picks[-1]
                    lon1 = lon0 + BAND - 0.5
                    yield _path(year, month, lon0, t0), _url(t0, t1, -77.5, 77.5, lon0, lon1, STRIDE), h
                    for name, (la0, la1, lo0, lo1) in DETAIL.items():
                        a, b = max(lo0, lon0), min(lo1, lon1)
                        if a <= b:
                            yield _path(year, month, lon0, t0, name), _url(t0, t1, la0, la1, a, b, 1), h


def _expected(month, region):
    """Paths a DETAIL box needs for one calendar month."""
    return [p for p, _, _ in _jobs() if p.name[5:].startswith(f"{month:02d}_{region}_")]


def fetch(workers=4):
    jobs = [j for j in _jobs() if not j[0].exists()]
    print(f"WW3: {len(jobs)} files to download", flush=True)

    def one(job):
        path, url, h = job
        _download(url, path)
        with xr.open_dataset(path) as ds:
            hours = set(ds.time.dt.hour.values.tolist())
        if hours != {h}:
            path.unlink()
            raise RuntimeError(f"{path.name}: expected hour {h}, got {sorted(hours)}")
        print(f"  {path.name}", flush=True)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, jobs))


def _open(path):
    ds = xr.open_dataset(path).squeeze("depth", drop=True)
    ds = ds.assign_coords(day=ds.time.dt.floor("D")).swap_dims(time="day").drop_vars("time")
    return ds[list(VARS)].load()


def _combine(bands):
    ds = xr.concat(bands, "longitude", join="outer")
    ds = ds.assign_coords(longitude=((ds.longitude + 180) % 360) - 180).sortby("longitude")
    return ds.rename(latitude="lat", longitude="lon")


def load(month, region=None):
    """Daily morning samples for one calendar month, all its seasons: dims (day, lat, lon).

    The global 1 degree grid by default, or a DETAIL box at 0.5 degrees. Longitudes are
    converted to -180..180. Bands can differ in which days exist (outages fall on
    different UTC hours), so missing band-days are NaN. Returns None for a DETAIL box
    whose files for the month are not all downloaded yet.
    """
    tag = f"_{region}" if region else ""
    if region and not all(p.exists() for p in _expected(month, region)):
        return None
    bands = []
    for lon0 in range(0, 360, BAND):
        parts = [_open(p) for year in config.seasons(month)
                 for p in sorted(OUT.glob(f"{year}-{month:02d}{tag}_lon{lon0:03d}_d*.nc"))]
        if parts:
            bands.append(xr.concat(parts, "day"))
    return _combine(bands)


def days(month):
    """Every UTC date in the month that has at least one model hour."""
    return sorted({t.date() for y in config.seasons(month) for t in available_hours(y, month)})
