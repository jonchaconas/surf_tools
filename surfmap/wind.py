"""Morning 10 m wind at each named break and stretch, from NOAA GFS (PacIOOS ERDDAP).

GFS here is 0.5 degree and 3-hourly, and its archive starts 2022-12-01, so the wind
covers fewer seasons than the waves. Each point is sampled OFFSHORE_DEG out to sea along
the direction it faces, at the 3-hourly step nearest config.LOCAL_HOUR solar time.
"""

import csv
import io
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import date, datetime, timedelta
from functools import lru_cache

import numpy as np
import requests

from . import config
from .spots import SPOTS
from .stretches import STRETCHES

ERDDAP = "https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ncep_global"
START = datetime(2022, 12, 1, 12)  # first stamp in the archive
EVERY_DAYS = 3                     # every third morning: ~475 per point, and 3x quicker to fetch
OFFSHORE_DEG = 0.3
KT = 1.943844                      # knots per m/s
OUT = config.DATA / "wind"


def _cell(lat, lon, facing):
    """The GFS grid point OFFSHORE_DEG out to sea from a coast point."""
    b = np.radians(facing)
    la = lat + OFFSHORE_DEG * np.cos(b)
    lo = lon + OFFSHORE_DEG * np.sin(b) / np.cos(np.radians(lat))
    return round(la * 2) / 2, (round(lo * 2) / 2) % 360


def _hour(lon):
    """The 3-hourly UTC step nearest config.LOCAL_HOUR solar time."""
    return int(round((config.LOCAL_HOUR - lon / 15) / 3) * 3) % 24


def _key(lat, lon, facing):
    """(grid lat, grid lon, UTC hour): breaks sharing a GFS cell and hour share one download."""
    return (*_cell(lat, lon, facing), _hour(lon))


def _path(key):
    la, lo, h = key
    return OUT / f"{la:+.1f}_{lo:05.1f}_{h:02d}z.csv"


def _url(key):
    la, lo, h = key
    # This dataset's latitude axis is mirrored: asking for 39.5 returns 39.5S. Checked against
    # ERA5 (Open-Meteo archive) on matching mornings, so request the opposite sign.
    la = -la
    t0 = START.replace(hour=h) + timedelta(days=0 if h >= START.hour else 1)
    last = datetime.combine(date.today() - timedelta(days=1), datetime.min.time()).replace(hour=h)
    t1 = t0 + timedelta(days=(last - t0).days // EVERY_DAYS * EVERY_DAYS)
    sel = f"[({t0:%Y-%m-%dT%H}:00:00Z):{8 * EVERY_DAYS}:({t1:%Y-%m-%dT%H}:00:00Z)][({la})][({lo})]"
    q = ",".join(f"{v}{sel}" for v in ("ugrd10m", "vgrd10m"))
    return f"{ERDDAP}.csv?" + q.replace("[", "%5B").replace("]", "%5D")


@lru_cache
def _keys():
    """Point name -> download key, for every named break, stretch and regional coast."""
    from .regions import REGIONS

    coasts = [c for spec in REGIONS.values() for c in spec["coasts"]]
    return {name: _key(lat, lon, facing) for name, _, lat, lon, facing, *_ in SPOTS + STRETCHES + coasts}


def fetch(workers=4):
    OUT.mkdir(parents=True, exist_ok=True)
    todo = sorted({k for k in _keys().values() if not _path(k).exists()})
    print(f"wind: {len(todo)} grid points to download", flush=True)

    def one(key):
        for attempt in range(4):
            try:
                r = requests.get(_url(key), timeout=600)
                r.raise_for_status()
                break
            except requests.RequestException:
                if attempt == 3:
                    raise
                time.sleep(15 * (attempt + 1))
        _path(key).write_text(r.text)
        print(f"  {_path(key).name}", flush=True)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(one, todo))


@lru_cache(maxsize=None)
def mornings(name):
    """(dates, speed in knots, direction the wind blows from) for one point's mornings."""
    rows = list(csv.reader(io.StringIO(_path(_keys()[name]).read_text())))[2:]
    t = np.array([r[0][:10] for r in rows], dtype="datetime64[D]")
    u, v = (np.array([float(r[i]) for r in rows]) for i in (3, 4))
    return t, np.hypot(u, v) * KT, np.degrees(np.arctan2(-u, -v)) % 360


def typical(name, months):
    """(median knots, most common of 16 directions) over the given calendar months, or None."""
    if name not in _keys() or not _path(_keys()[name]).exists():
        return None
    t, kt, frm = mornings(name)
    sel = np.isin(t.astype("datetime64[M]").astype(int) % 12 + 1, months) & np.isfinite(kt)
    if not sel.any():
        return None
    sector = np.round(frm[sel] / 22.5).astype(int) % 16
    return float(np.median(kt[sel])), float(np.bincount(sector, minlength=16).argmax() * 22.5)
