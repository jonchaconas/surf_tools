"""Morning wind along every coast, for the map's coastal colours (the marked spots have their own).

ERA5 10 m wind (CDS point time series) at each coast point (export._coast_points: one per ~200 km
of shoreline where nothing else is mapped) and each OpenStreetMap surf spot, sampled at every
config.HOURS local solar time, 2017 on. Together with the marked spots that's a wind series within
~100 km of any shore. Slow: the CDS takes about two requests at a time, so ~2,700 points run
overnight; re-running resumes.
"""

import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, timedelta

import numpy as np
import xarray as xr

from . import config

OUT = config.DATA / "coastwind"
WEB = config.ROOT / "web" / "data"


def _path(lat, lon):
    return OUT / f"{lat:+.3f}_{lon:+.3f}.npz"


def points():
    """[(lat, lon)] of every coast point and OpenStreetMap spot the atlas shows."""
    gaps = json.loads((WEB / "coastpts.json").read_text())
    osm = json.loads((WEB / "osm.json").read_text())
    return sorted({(round(g[2], 3), round(g[3], 3)) for g in gaps} | {(round(o[1], 3), round(o[2], 3)) for o in osm})


def _fetch(lat, lon):
    """ERA5 u/v at one point, kept only at the config.HOURS local times: npz of day, hour, kt, from."""
    import tempfile
    import zipfile

    import cdsapi

    from .points import KT, _CDS, _local_mornings

    path = _path(lat, lon)
    if path.exists():
        return
    with tempfile.TemporaryDirectory() as tmp:
        dl = f"{tmp}/w.zip"
        request = {"variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
                   "location": {"latitude": lat, "longitude": lon},
                   "date": [f"2016-12-31/{date.today() - timedelta(days=6):%Y-%m-%d}"], "data_format": "netcdf"}
        for attempt in range(6):
            try:
                with _CDS:
                    cdsapi.Client(quiet=True, progress=False).retrieve("reanalysis-era5-single-levels-timeseries", request).download(dl)
                break
            except Exception as e:
                if attempt == 5 or "queued requests" not in str(e):
                    raise
                import time
                time.sleep(60)
        if zipfile.is_zipfile(dl):
            with zipfile.ZipFile(dl) as z:
                z.extractall(tmp)
                nc = next(f"{tmp}/{n}" for n in z.namelist() if n.endswith(".nc"))
        else:
            nc = dl
        with xr.open_dataset(nc) as d:
            t = d["valid_time"].values if "valid_time" in d else d["time"].values
            u, v = d["u10"].values.astype(float), d["v10"].values.astype(float)
    days, hours, kts, froms = [], [], [], []
    for h in config.HOURS:
        idx, dd = _local_mornings(t, lon, 1, h)
        days.append(dd); hours.append(np.full(len(dd), h, np.int8))
        kts.append(np.hypot(u[idx], v[idx]) * KT); froms.append(np.degrees(np.arctan2(-u[idx], -v[idx])) % 360)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, day=np.concatenate(days).astype("datetime64[D]").astype(np.int32), hour=np.concatenate(hours),
                        kt=np.concatenate(kts).astype(np.float32), wind_from=np.concatenate(froms).astype(np.float32))


def fetch(workers=2):
    todo = [p for p in points() if not _path(*p).exists()]
    print(f"coast wind: {len(todo)} points to fetch", flush=True)
    done = failed = 0
    with ThreadPoolExecutor(workers) as pool:
        futures = {pool.submit(_fetch, *p): p for p in todo}
        for fut in as_completed(futures):
            try:
                fut.result(); done += 1
            except Exception as e:
                failed += 1
                print(f"  {futures[fut]} failed: {e}", flush=True)
            if (done + failed) % 50 == 0:
                print(f"  {done + failed}/{len(todo)} ({failed} failed)", flush=True)
    print(f"done: {done}, failed: {failed}", flush=True)
