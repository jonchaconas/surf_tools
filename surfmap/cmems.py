"""Fill Open-Meteo gaps from Copernicus Marine, which serves the same Meteo-France MFWAM model.

Open-Meteo's free tier caps requests per day. Copernicus has MFWAM's global analysis
(ANFC, 1/12 degree, 3-hourly, from November 2022) and the WAVERYS reanalysis (REAN,
0.2 degree, 3-hourly) for earlier months. This writes any inland.py or regions.py file
that is still missing, in exactly the layout those modules read, so nothing downstream
changes. Copernicus has no swell peak period; readers estimate it from the mean swell
period, as for the Great Lakes. Needs a Copernicus Marine login (`copernicusmarine login`).
"""

import json
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta

import copernicusmarine as cm
import numpy as np
import xarray as xr

from . import config, inland, regions

ANFC = "cmems_mod_glo_wav_anfc_0.083deg_PT3H-i"
REAN = "cmems_mod_glo_wav_my_0.2deg_PT3H-i"
ANFC_START = datetime(2022, 11, 1)
SPAN = (datetime(2021, 10, 1), datetime(2026, 11, 1))
VARS = {  # Open-Meteo name -> Copernicus variable
    "swell_wave_height": "VHM0_SW1", "swell_wave_period": "VTM01_SW1", "swell_wave_direction": "VMDR_SW1",
    "wind_wave_height": "VHM0_WW", "wave_height": "VHM0", "wave_period": "VTM10", "wave_direction": "VMDR",
}
PAD = 0.3  # degrees searched around each point for the nearest cell with data
OUT = config.DATA / "cmems"
_WRITE = threading.Lock()  # HDF5 is not thread-safe: one netCDF write at a time


def _missing():
    """(path, [(lat, lon), ...], date window) for every inland/regional file not yet on disk."""
    jobs = []
    for month in config.MONTHS:
        for year in config.seasons(month):
            d0, d1 = inland._window(year, month)
            for name in inland.BASINS:
                lats, lons = inland.points(name)
                for k, path in enumerate(inland._files(name, year, month)):
                    if not path.exists():
                        pts = list(zip(lats[k * inland.CHUNK:(k + 1) * inland.CHUNK],
                                       lons[k * inland.CHUNK:(k + 1) * inland.CHUNK]))
                        jobs.append((path, [(float(a), float(b)) for a, b in pts], d0, d1))
            for region, spec in regions.REGIONS.items():
                path = regions._path(region, year, month)
                if not path.exists():
                    jobs.append((path, [(c[2], c[3]) for c in spec["coasts"]], d0, d1))
    return jobs


def _cache(lat, lon):
    return OUT / f"{lat:+.3f}_{lon:+.3f}.nc"


def _nearest_series(dataset, lat, lon, t0, t1):
    """3-hourly series at the cell nearest (lat, lon) that has waves, from one dataset."""
    ds = cm.open_dataset(dataset_id=dataset, variables=list(VARS.values()),
                         minimum_longitude=lon - PAD, maximum_longitude=lon + PAD,
                         minimum_latitude=lat - PAD, maximum_latitude=lat + PAD,
                         start_datetime=t0.isoformat(), end_datetime=t1.isoformat())
    sea = np.isfinite(ds["VHM0"].isel(time=slice(0, 16)).values).any(axis=0)
    la, lo = np.meshgrid(ds.latitude.values, ds.longitude.values, indexing="ij")
    km = np.hypot(la - lat, (lo - lon) * np.cos(np.radians(lat)))
    i, j = np.unravel_index(np.argmin(np.where(sea, km, np.inf)), km.shape)
    return ds.isel(latitude=i, longitude=j).load()


def _fetch_point(point):
    lat, lon = point
    path = _cache(lat, lon)
    if path.exists():
        return
    parts = [_nearest_series(REAN, lat, lon, SPAN[0], ANFC_START - timedelta(hours=3)),
             _nearest_series(ANFC, lat, lon, ANFC_START, SPAN[1])]
    ds = xr.concat([p.drop_vars(["latitude", "longitude"]).assign_attrs(cell=(float(p.latitude), float(p.longitude)))
                    for p in parts], "time")
    tmp = path.with_suffix(".part")
    with _WRITE:
        ds.to_netcdf(tmp)
    tmp.rename(path)
    print(f"  {path.name}", flush=True)


def _rows(points, d0, d1):
    """Open-Meteo-shaped rows: hourly times, each hour taking the nearest 3-hourly value."""
    hours = np.arange(np.datetime64(d0, "h"), np.datetime64(d1, "h") + 24, np.timedelta64(1, "h"))
    rows = []
    for lat, lon in points:
        with xr.open_dataset(_cache(lat, lon)) as ds:
            t = ds.time.values.astype("datetime64[h]")
            idx = np.clip(np.searchsorted(t, hours), 1, len(t) - 1)
            idx -= (hours - t[idx - 1]) < (t[idx] - hours)
            ok = np.abs(t[idx] - hours) <= np.timedelta64(2, "h")
            hourly = {"time": [str(h)[:13] + ":00" for h in hours]}
            for om, cv in VARS.items():
                v = ds[cv].values[idx].astype(float)
                hourly[om] = [None if (not k or np.isnan(x)) else round(float(x), 2) for x, k in zip(v, ok)]
            hourly["swell_wave_peak_period"] = [None] * len(hours)
        rows.append({"requested": [lat, lon], "hourly": hourly, "source": "copernicus"})
    return rows


def fill(workers=4):
    jobs = _missing()
    points = sorted({p for _, pts, _, _ in jobs for p in pts})
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Copernicus: {len(jobs)} files missing, {len(points)} points to fetch", flush=True)
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(_fetch_point, points))
    for path, pts, d0, d1 in jobs:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(_rows(pts, d0, d1)))
    print(f"wrote {len(jobs)} files from Copernicus", flush=True)
