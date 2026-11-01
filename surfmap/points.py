"""Morning sea state and wind at each named break, stretch and regional coast.

Waves come from Copernicus Marine's MFWAM at the sea cell nearest each point, a little out
to sea along the way it faces: the 1/12 degree analysis from November 2022, and the 0.2 degree
WAVERYS reanalysis before that. Both swell partitions are kept, so a small groundswell running
under a local windswell (Malibu's summer south swells, say) is still seen. Wind is ERA5's 10 m
wind at the grid point just offshore, from the Climate Data Store's point time-series dataset.
Each morning is the step nearest config.LOCAL_HOUR solar time, dated by the local day.

Needs a Copernicus Marine login (`copernicusmarine login`) and a CDS key in ~/.cdsapirc.
"""

import threading
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, timedelta
from functools import lru_cache

import numpy as np
import xarray as xr

from . import config, regions, tides
from .spots import SPOTS
from .stretches import STRETCHES

OUT = config.DATA / "points"
WAVES = ["VHM0", "VTM10", "VMDR", "VHM0_SW1", "VTM01_SW1", "VMDR_SW1",
         "VHM0_SW2", "VTM01_SW2", "VMDR_SW2", "VHM0_WW"]
SEAWARD_DEG = 0.05   # pick the wave cell nearest a spot this far out along its facing
WIND_DEG = 0.15      # and the ERA5 point this far out, so it is over water
KT = 1.943844        # knots per m/s
_WRITE = threading.Lock()  # HDF5 is not thread-safe
_CDS = threading.BoundedSemaphore(2)  # the CDS rejects bigger queues per dataset ("request sequentially")


def all_points():
    """[(kind, name, region, lat, lon, facing)] for every named break, stretch and regional coast."""
    out = [("break", n, r, la, lo, f) for n, r, la, lo, f in SPOTS]
    out += [("stretch", n, r, la, lo, f) for n, r, la, lo, f, _ in STRETCHES]
    out += [("coast", c[0], region, c[2], c[3], c[4]) for region, spec in regions.REGIONS.items() for c in spec["coasts"]]
    return out


def _out(lat, lon, facing, deg):
    b = np.radians(facing)
    return lat + deg * np.cos(b), lon + deg * np.sin(b) / np.cos(np.radians(lat))


def _stem(lat, lon):
    return f"{lat:+.3f}_{lon:+.3f}"


def _waves_path(lat, lon, early=False):
    return OUT / f"{_stem(lat, lon)}_waves{'_early' if early else ''}.nc"


def _wind_path(lat, lon, early=False):
    return OUT / f"{_stem(lat, lon)}_wind{'_early' if early else ''}.nc"


def _nearest(dataset, lat, lon, facing, t0, t1, pad=0.3):
    import copernicusmarine as cm

    ds = cm.open_dataset(dataset_id=dataset, variables=WAVES,
                         minimum_longitude=lon - pad, maximum_longitude=lon + pad,
                         minimum_latitude=lat - pad, maximum_latitude=lat + pad,
                         start_datetime=t0.isoformat(), end_datetime=t1.isoformat())
    sea = np.isfinite(ds["VHM0"].isel(time=slice(0, 16)).values).any(axis=0)
    tla, tlo = _out(lat, lon, facing, SEAWARD_DEG)
    la, lo = np.meshgrid(ds.latitude.values, ds.longitude.values, indexing="ij")
    km = np.hypot(la - tla, (lo - tlo) * np.cos(np.radians(lat)))
    i, j = np.unravel_index(np.argmin(np.where(sea, km, np.inf)), km.shape)
    if not sea[i, j]:
        raise ValueError(f"no sea cell within {pad} degrees of {lat}, {lon}")
    return ds.isel(latitude=i, longitude=j).load()


def _fetch_waves(lat, lon, facing, early=False):
    """Copernicus waves: from cmems.SPAN on, or with early=True the reanalysis years before it."""
    from .cmems import ANFC, ANFC_START, EARLY, REAN, SPAN  # imports the Copernicus client: only when fetching

    path = _waves_path(lat, lon, early)
    if path.exists():
        return
    parts = [_nearest(REAN, lat, lon, facing, *EARLY)] if early else \
        [_nearest(REAN, lat, lon, facing, SPAN[0], ANFC_START - timedelta(hours=3)),
         _nearest(ANFC, lat, lon, facing, ANFC_START, SPAN[1])]
    ds = xr.concat([p.drop_vars(["latitude", "longitude"]) for p in parts], "time")
    ds.attrs["cell"] = (float(parts[-1].latitude), float(parts[-1].longitude))  # the 1/12 degree cell
    tmp = path.with_suffix(".part")
    with _WRITE:
        ds.to_netcdf(tmp)
    tmp.rename(path)


def _fetch_wind(lat, lon, facing, early=False):
    import cdsapi

    from .cmems import EARLY, SPAN

    path = _wind_path(lat, lon, early)
    if path.exists():
        return
    wla, wlo = _out(lat, lon, facing, WIND_DEG)
    first, last = (EARLY[0], EARLY[1]) if early else (SPAN[0], date.today() - timedelta(days=6))  # ERA5 is ~5 days behind
    tmp = path.with_suffix(".zip")
    request = {
        "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
        "location": {"latitude": round(float(wla), 3), "longitude": round(float(((wlo + 180) % 360) - 180), 3)},
        "date": [f"{first:%Y-%m-%d}/{last:%Y-%m-%d}"],
        "data_format": "netcdf",
    }
    for attempt in range(6):
        try:
            with _CDS:
                cdsapi.Client(quiet=True, progress=False).retrieve("reanalysis-era5-single-levels-timeseries", request).download(str(tmp))
            break
        except Exception as e:  # a full queue clears in a minute or two
            if attempt == 5 or "queued requests" not in str(e):
                raise
            time.sleep(60)
    if zipfile.is_zipfile(tmp):
        with zipfile.ZipFile(tmp) as z:
            data = z.read(next(n for n in z.namelist() if n.endswith(".nc")))
        path.with_suffix(".part").write_bytes(data)
        tmp.unlink()
    else:
        tmp.rename(path.with_suffix(".part"))
    path.with_suffix(".part").rename(path)


def fetch(workers=4, early=False):
    """Download waves (Copernicus) and wind (ERA5) for every point not yet cached; early=True
    fetches the years before cmems.SPAN, back to the WW3 archive's start."""
    OUT.mkdir(parents=True, exist_ok=True)
    pts = {(la, lo): f for _, _, _, la, lo, f in all_points()}
    jobs = [(fn, la, lo, f) for (la, lo), f in pts.items()
            for fn, path in ((_fetch_waves, _waves_path), (_fetch_wind, _wind_path)) if not path(la, lo, early).exists()]
    print(f"points: {len(pts)}, downloads to do: {len(jobs)}", flush=True)
    done = failed = 0
    with ThreadPoolExecutor(workers) as pool:
        futures = {pool.submit(fn, la, lo, f, early): (fn.__name__, la, lo) for fn, la, lo, f in jobs}
        for fut in as_completed(futures):
            what, la, lo = futures[fut]
            try:
                fut.result()
                done += 1
            except Exception as e:  # keep going; a re-run picks up what failed
                failed += 1
                print(f"  {what} {la},{lo} failed: {e}", flush=True)
            if (done + failed) % 20 == 0:
                print(f"  {done + failed}/{len(jobs)} ({failed} failed)", flush=True)
    print(f"done: {done}, failed: {failed}", flush=True)


def _local_mornings(t_utc, lon, step_h):
    """Indices of the step nearest config.LOCAL_HOUR solar time on each local day, and those days."""
    local = t_utc + np.timedelta64(int(round(lon / 15 * 60)), "m")
    hours = (local - local.astype("datetime64[D]")) / np.timedelta64(1, "h")
    keep = np.flatnonzero(np.abs(hours - config.LOCAL_HOUR) <= step_h / 2)
    days = local[keep].astype("datetime64[D]")
    first = np.unique(days, return_index=True)[1]
    return keep[first], days[first]


def _open(paths):
    """One dataset over the files that exist, in time order (the early years' file first)."""
    paths = [p for p in paths if p.exists()]
    if len(paths) == 1:
        return xr.open_dataset(paths[0])
    dims = [d for d in ("time", "valid_time") if d in xr.open_dataset(paths[0]).dims]
    return xr.concat([xr.open_dataset(p) for p in paths], dims[0], data_vars="minimal", coords="minimal", compat="override")


@lru_cache(maxsize=None)
def mornings(lat, lon):
    """Morning arrays for one point, or None until both downloads are cached:
    day, hs1/tp1/dir1 and hs2/tp2/dir2 (swell partitions; tp estimated from the mean period),
    hs/tm/dir (all waves), ww (wind-sea height), kt and wind_from (ERA5 10 m wind), tide (m, EOT20), cell, source.
    In California the waves come from CDIP's nearshore model (cdip.py) instead."""
    wp, vp = _waves_path(lat, lon), _wind_path(lat, lon)
    if not (wp.exists() and vp.exists()):
        return None
    from .inland import TP_RATIO

    with _open([_waves_path(lat, lon, True), wp]) as w:
        idx, days = _local_mornings(w.time.values, lon, 3)
        g = lambda v: w[v].values[idx].astype(float)
        tide = tides.height(lat, lon, w.time.values[idx])
        out = dict(day=days, cell=tuple(w.attrs["cell"]), tide=np.full(len(days), np.nan) if tide is None else tide,
                   hs1=g("VHM0_SW1"), tp1=g("VTM01_SW1") * TP_RATIO, dir1=g("VMDR_SW1"),
                   hs2=g("VHM0_SW2"), tp2=g("VTM01_SW2") * TP_RATIO, dir2=g("VMDR_SW2"),
                   hs=g("VHM0"), tm=g("VTM10"), dir=g("VMDR"), ww=g("VHM0_WW"))
    with _open([_wind_path(lat, lon, True), vp]) as v:
        t = v["valid_time"].values if "valid_time" in v else v["time"].values
        widx, wdays = _local_mornings(t, lon, 1)
        u, vv = v["u10"].values[widx].astype(float), v["v10"].values[widx].astype(float)
    pos = np.searchsorted(wdays, days)
    hit = (pos < len(wdays)) & (wdays[np.minimum(pos, len(wdays) - 1)] == days)
    out["kt"] = np.where(hit, np.hypot(u, vv)[np.minimum(pos, len(wdays) - 1)] * KT, np.nan)
    out["wind_from"] = np.where(hit, (np.degrees(np.arctan2(-u, -vv)) % 360)[np.minimum(pos, len(wdays) - 1)], np.nan)
    out["source"] = "MFWAM ~8 km cell"

    # Bend the offshore swell trains in to the break over the seafloor (refract.py), losses only. Off
    # for now (config.REFRACT): it matches CDIP in California, but elsewhere ETOPO's seafloor is too
    # coarse for reefs, points and canyons and it all but erases Puerto Escondido, the Mentawais,
    # Chicama and Jeffreys Bay.
    from . import refract

    bend = refract.transfer(lat, lon, max_k=1.0) if config.REFRACT else None
    if bend is not None:
        for n in ("1", "2"):  # heights only: windows and exposure are judged on the offshore direction
            out["hs" + n] = bend(out["hs" + n], out["tp" + n], out["dir" + n])[0]
        out["source"] = "MFWAM ~8 km cell, refracted to the break"

    # California: CDIP's nearshore model replaces the offshore waves wherever it has the morning.
    from . import cdip

    near = cdip.mornings(lat, lon)
    if near is not None:
        pos = np.minimum(np.searchsorted(near["day"], days), len(near["day"]) - 1)
        hit = (near["day"][pos] == days) & np.isfinite(near["hs"][pos])
        for k in ("hs1", "tp1", "dir1", "hs2", "tp2", "dir2", "hs", "dir"):
            out[k] = np.where(hit, near[k][pos], out[k])
        out["tm"] = np.where(hit, near["tp"][pos] / TP_RATIO, out["tm"])  # read back as tm * TP_RATIO
        out["cell"] = near["cell"]
        out["source"] = f"CDIP MOP site {near['site']}, {near['depth']:.0f} m deep"
    return out


@lru_cache(maxsize=None)
def _days(month):
    from . import ww3

    return np.array(ww3.days(month), dtype="datetime64[D]")


def season(m, months):
    """The mornings of m in the calendar months' seasons that WW3 also covers (so every source
    counts the same days) and that have both waves and wind."""
    sel = np.isin(m["day"], np.concatenate([_days(mo) for mo in months])) & np.isfinite(m["hs"]) & np.isfinite(m["kt"])
    return {k: v[sel] if isinstance(v, np.ndarray) else v for k, v in m.items()}


def typical_wind(s):
    """(median knots, most common of 16 directions) over the mornings in s, or None."""
    if not len(s["kt"]):
        return None
    sector = np.round(s["wind_from"] / 22.5).astype(int) % 16
    return float(np.median(s["kt"])), float(np.bincount(sector, minlength=16).argmax() * 22.5)


def summary(lat, lon, facing, months, relaxed=False, window=None):
    """One point over the calendar months: % good, mornings, the median good morning (hs, hb, tp,
    dir) and the typical wind, as analysis._merge; None until it is downloaded. relaxed applies
    the enclosed-sea standard to total waves (no wind test), as the map insets do."""
    from . import rating
    from .analysis import _circmean, _relaxed
    from .inland import TP_RATIO

    m = mornings(lat, lon)
    if m is None:
        return None
    s = season(m, months)
    if relaxed:
        good = _relaxed(s["hs"], s["tm"], TP_RATIO) & rating.exposed(s["dir"], facing)
        hs, tp, d = s["hs"], s["tm"] * TP_RATIO, s["dir"]
    else:
        good = rating.break_good(s, facing, window=window)
        hs, tp, d = rating.good_swell(s, facing, window)
    hs, tp, d = hs[good], tp[good], d[good]
    med = lambda a: float(np.median(a)) if len(a) else np.nan
    return dict(pct=100 * good.sum() / max(len(good), 1), days=len(good), hs=med(hs),
                hb=med(rating.breaker_height(hs, tp)), tp=med(tp), dir=_circmean(d), wind=typical_wind(s))
