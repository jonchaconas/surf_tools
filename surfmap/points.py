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

# MFWAM gives mean periods; surfers (and buoys) quote the peak. Calibrated on NDBC dominant periods
# at the buoys within 60 km of a point, 2017 on: a swell train's peak runs 1.28x its Tm01 at 7-8 s
# but only 1.08x past 15 s (a flat 1.29 made Pascuales' small forerunners 25-30 s), and the whole
# sea's 1.13x its Tm-1,0 in the ocean but 1.33x in the Great Lakes, where the model's periods run
# short. Enclosed seas, with no buoys here, take the lakes' figure.
SWELL_TP = ([7.5, 10, 12, 14, 16.5], [1.28, 1.25, 1.20, 1.13, 1.08])  # (Tm01 s, peak / Tm01)
SEA_TP, LAKE_SEA_TP = 1.13, 1.33  # peak / Tm-1,0
MAX_TP = 22.0      # longer than buoys report; a model forerunner
MIN_TRAIN_M = 0.2  # smaller swell partitions are model remnants, not a train to surf


def swell_tp(tm01):
    """A swell partition's peak period from its mean period (SWELL_TP), capped at MAX_TP."""
    tm01 = np.asarray(tm01, float)
    return np.minimum(tm01 * np.interp(tm01, *SWELL_TP), MAX_TP)


def sea_tp(tm10, lake):
    """The whole sea's peak period from its Tm-1,0, in the open ocean or (lake) a lake or enclosed sea."""
    return np.minimum(np.asarray(tm10, float) * (LAKE_SEA_TP if lake else SEA_TP), MAX_TP)


# Wave heights corrected to the buoys (buoys.quantile_maps()): each Great Lake by all its well-matched
# buoys together, and an ocean break by its own nearest well-matched buoy, if it has one. Quantile
# mapping keeps the model's day-to-day order but gives its heights the buoys' distribution, which
# matters most in the lakes, where MFWAM under-calls 4 ft days two to three times over.
QMAP = config.ROOT / "surfmap" / "buoy_qmap.json"  # small, and committed: the daily forecast job needs it
LAKES = {  # (lat0, lat1, lon0, lon1); the first box that holds a point names its lake
    "Lake Superior": (46.2, 49.1, -92.3, -84.3),
    "Lake Michigan": (41.5, 46.1, -88.1, -84.8),
    "Lake Huron": (43.0, 46.4, -84.8, -79.6),
    "Lake Erie": (41.3, 42.95, -83.6, -78.8),
    "Lake Ontario": (43.1, 44.3, -79.9, -76.0),
}


def lake_of(lat, lon):
    return next((name for name, (a, b, c, d) in LAKES.items() if a <= lat <= b and c <= lon <= d), None)


@lru_cache(maxsize=1)
def _qmaps():
    import json

    return json.loads(QMAP.read_text()) if QMAP.exists() else {}


def correction(lat, lon):
    """(model quantiles, buoy quantiles, what they come from) for this point, or None."""
    q, lake = _qmaps(), lake_of(lat, lon)
    e = q.get("lakes", {}).get(lake) if lake else q.get("points", {}).get(f"{lat:.3f},{lon:.3f}")
    return None if not e else (np.array(e["model"]), np.array(e["buoy"]), f"{lake} buoys" if lake else f"buoy {e['station']}")


def corrected(hs, qm):
    """Heights through a quantile map: interpolated within it, scaled by its end ratios beyond."""
    m, b = qm[0] + np.arange(len(qm[0])) * 1e-6, qm[1]  # strictly increasing for np.interp
    hs = np.asarray(hs, float)
    out = np.interp(hs, m, b)
    return np.where(hs > m[-1], hs * b[-1] / m[-1], np.where(hs < m[0], hs * b[0] / max(m[0], 1e-6), out))


def trains(hs1, tm1, dir1, hs2, tm2, dir2):
    """The two swell partitions with peak periods, partitions under MIN_TRAIN_M dropped (NaN)."""
    out = {}
    for n, (hs, tm, d) in (("1", (hs1, tm1, dir1)), ("2", (hs2, tm2, dir2))):
        keep = np.asarray(hs, float) >= MIN_TRAIN_M
        out["hs" + n], out["tp" + n], out["dir" + n] = (np.where(keep, v, np.nan) for v in (hs, swell_tp(tm), d))
    return out


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


def _local_mornings(t_utc, lon, step_h, hour=config.LOCAL_HOUR):
    """Indices of the step nearest `hour` local solar time on each local day, and those days."""
    local = t_utc + np.timedelta64(int(round(lon / 15 * 60)), "m")
    hours = (local - local.astype("datetime64[D]")) / np.timedelta64(1, "h")
    keep = np.flatnonzero(np.abs(hours - hour) <= step_h / 2)
    days = local[keep].astype("datetime64[D]")
    first = np.unique(days, return_index=True)[1]
    return keep[first], days[first]


def _open(paths):
    """One dataset over the files that exist, in time order (the early years' file first)."""
    paths = [p for p in paths if p.exists()]
    if len(paths) == 1:
        return xr.open_dataset(paths[0])
    dims = [d for d in ("time", "valid_time") if d in xr.open_dataset(paths[0]).dims]
    ds = xr.concat([xr.open_dataset(p) for p in paths], dims[0], data_vars="minimal", coords="minimal", compat="override")
    ds.attrs = xr.open_dataset(paths[-1]).attrs  # the current grid's cell, which the forecast reads, not the early 0.2° one
    return ds


@lru_cache(maxsize=None)
def mornings(lat, lon, hour=config.LOCAL_HOUR, raw=False):
    """Arrays for one point, a sample a day nearest `hour` local solar time (the morning by default),
    or None until both downloads are cached:
    day, hs1/tp1/dir1 and hs2/tp2/dir2 (swell partitions, peak period from the mean: swell_tp),
    hs/tm/tp/dir (all waves; tm the model's Tm-1,0, tp the peak: sea_tp), ww (wind-sea height), kt and wind_from (ERA5 10 m wind), tide (m, EOT20), cell, source.
    In California the waves come from CDIP's nearshore model (cdip.py) instead. Heights are corrected
    to nearby buoys where there are good ones (correction()), unless raw."""
    wp, vp = _waves_path(lat, lon), _wind_path(lat, lon)
    if not (wp.exists() and vp.exists()):
        return None
    from .analysis import _basin_of

    with _open([_waves_path(lat, lon, True), wp]) as w:
        idx, days = _local_mornings(w.time.values, lon, 3, hour)
        g = lambda v: w[v].values[idx].astype(float)
        tide = tides.height(lat, lon, w.time.values[idx])
        out = dict(day=days, cell=tuple(w.attrs["cell"]), tide=np.full(len(days), np.nan) if tide is None else tide,
                   **trains(g("VHM0_SW1"), g("VTM01_SW1"), g("VMDR_SW1"), g("VHM0_SW2"), g("VTM01_SW2"), g("VMDR_SW2")),
                   hs=g("VHM0"), tm=g("VTM10"), tp=sea_tp(g("VTM10"), _basin_of(lat, lon)), dir=g("VMDR"), ww=g("VHM0_WW"))
    with _open([_wind_path(lat, lon, True), vp]) as v:
        t = v["valid_time"].values if "valid_time" in v else v["time"].values
        widx, wdays = _local_mornings(t, lon, 1, hour)
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

    near = cdip.mornings(lat, lon, hour)
    if near is not None:
        pos = np.minimum(np.searchsorted(near["day"], days), len(near["day"]) - 1)
        hit = (near["day"][pos] == days) & np.isfinite(near["hs"][pos])
        for k in ("hs1", "tp1", "dir1", "hs2", "tp2", "dir2", "hs", "tp", "dir"):
            out[k] = np.where(hit, near[k][pos], out[k])
        out["tm"] = np.where(hit, np.nan, out["tm"])  # CDIP gives the peak period only
        out["cell"] = near["cell"]
        out["source"] = f"CDIP MOP site {near['site']}, {near['depth']:.0f} m deep"

    qm = None if raw else correction(lat, lon)
    if qm is not None:  # every height by the same factor as the whole sea's
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(out["hs"] > 0, corrected(out["hs"], qm) / out["hs"], 1.0)
        for k in ("hs", "hs1", "hs2", "ww"):
            out[k] = out[k] * f
        out["source"] += f", heights corrected to {qm[2]}"
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

    m = mornings(lat, lon)
    if m is None:
        return None
    s = season(m, months)
    if relaxed:
        good = _relaxed(s["hs"], s["tp"], 1.0) & rating.exposed(s["dir"], facing)
        hs, tp, d = s["hs"], s["tp"], s["dir"]
    else:
        good = rating.break_good(s, facing, window=window)
        hs, tp, d = rating.good_swell(s, facing, window)
    hs, tp, d = hs[good], tp[good], d[good]
    med = lambda a: float(np.median(a)) if len(a) else np.nan
    return dict(pct=100 * good.sum() / max(len(good), 1), days=len(good), hs=med(hs),
                hb=med(rating.breaker_height(hs, tp)), tp=med(tp), dir=_circmean(d), wind=typical_wind(s))
