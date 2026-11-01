"""Monthly sea-surface temperature at any coast, for the wetsuit guide.

NOAA OISST v2.1 (1/4 degree), 1991-2020 long-term monthly means, one 46 MB file from NOAA PSL with
no login: https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
"""

from functools import lru_cache

import numpy as np
import requests
import xarray as xr

from . import config

URL = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.mon.ltm.1991-2020.nc"
PATH = config.DATA / "sst" / "sst.mon.ltm.1991-2020.nc"
SEARCH = 0.75  # degrees: how far to look for the nearest cell with water

# Wetsuit by the water temperature (deg C, at least): a common surf-shop guide.
SUITS = [(24, "boardshorts"), (22, "rash vest or shorty"), (19, "2 mm shorty or 3/2"), (16, "3/2 full"),
         (13, "4/3 full"), (10, "5/4 full, boots"), (-5, "6/5 hooded, boots, gloves")]


def fetch():
    if PATH.exists():
        return
    PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = PATH.with_suffix(".part")
    with requests.get(URL, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    tmp.rename(PATH)


@lru_cache
def _grid():
    with xr.open_dataset(PATH, decode_times=False) as ds:
        return ds.lat.values, ds.lon.values, ds.sst.values  # (12, lat, lon), deg C; NaN over land


def monthly(lat, lon):
    """Twelve monthly mean temperatures (deg C, Jan..Dec) at the nearest cell with water, or None."""
    if not PATH.exists():
        return None
    la, lo, sst = _grid()
    i = np.flatnonzero(np.abs(la - lat) <= SEARCH)
    j = np.flatnonzero(np.abs((lo - lon % 360 + 180) % 360 - 180) <= SEARCH / max(np.cos(np.radians(lat)), 0.2))
    if not len(i) or not len(j):
        return None
    box = sst[:, i][:, :, j]
    ok = np.isfinite(box).all(axis=0)
    if not ok.any():
        return None
    gi, gj = np.meshgrid(la[i], lo[j], indexing="ij")
    d = np.where(ok, np.hypot(gi - lat, ((gj - lon % 360 + 180) % 360 - 180) * np.cos(np.radians(lat))), np.inf)
    a, b = np.unravel_index(np.argmin(d), d.shape)
    return [round(float(v), 1) for v in box[:, a, b]]


def suit(celsius):
    return next(name for floor, name in SUITS if celsius >= floor)
