"""Nearshore waves at California breaks from CDIP's MOP model (Scripps Institution of Oceanography).

MOP (v1.1) carries offshore wave spectra in to 10 m depth every 100 m along the California coast,
so island shadows (the Channel Islands), refraction round headlands and shoaling are modelled,
which the ~8 km offshore models can't do. It is hourly from 2000 (hindcast, then nowcast) and
served over OPeNDAP with no login: https://cdip.ucsd.edu/MOP_v1.1/

Each break is matched to its nearest MOP site. Each morning's 1-D spectrum is split at the troughs
between its peaks into wave trains; the two biggest stand in for the offshore models' two swell
partitions, so the good+ rule is unchanged. points.mornings() swaps these in for the offshore
waves wherever a break has a site (ERA5 still supplies the wind).
"""

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

import numpy as np
import requests
import xarray as xr

from . import config

THREDDS = "https://thredds.cdip.ucsd.edu/thredds"
URL = THREDDS + "/dodsC/cdip/model/MOP_alongshore/{}_{}.nc"
BOX = (32.4, 42.1, -124.6, -117.0)  # California: breaks inside it are looked up
MAX_KM = 3.0                        # a break further than this from every site keeps the offshore model
COARSE = 25                         # every 25th site (~2.5 km) for the first pass of the lookup
OUT = config.DATA / "cdip"
SITES = OUT / "sites.json"


@lru_cache(maxsize=None)
def _catalog():
    """{prefix: [site IDs in order along the coast]}, from the THREDDS catalog (the numbers have gaps)."""
    path = OUT / "catalog.json"
    if not path.exists():
        xml = requests.get(THREDDS + "/catalog/cdip/model/MOP_alongshore/catalog.xml", timeout=120).text
        names = sorted(set(re.findall(r'name="([A-Z]+\d+)_hindcast\.nc"', xml)))
        OUT.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(names))
    out = {}
    for name in json.loads(path.read_text()):
        out.setdefault(name.rstrip("0123456789"), []).append(name)
    return out


def _meta(site):
    """(site, [lat, lon, shore normal]) from the OPeNDAP text endpoint, or (site, None)."""
    for attempt in range(4):
        try:
            r = requests.get(URL.format(site, "hindcast") + ".ascii?metaLatitude,metaLongitude,metaShoreNormal", timeout=30)
            if r.status_code == 404:
                return site, None
            r.raise_for_status()
            vals = [float(v) for v in re.findall(r"^meta\w+, ([-\d.eE]+)$", r.text, re.M)]
            return site, vals if len(vals) == 3 else None
        except requests.RequestException:
            time.sleep(2 * (attempt + 1))
    return site, None


def _sites(names):
    """{site: [lat, lon, shore normal]} for the given site IDs, cached in data/cdip/sites.json."""
    known = {k: v for k, v in (json.loads(SITES.read_text()) if SITES.exists() else {}).items() if v}
    todo = [n for n in names if n not in known]
    if todo:
        with ThreadPoolExecutor(6) as pool:
            known.update(pool.map(_meta, todo))
        OUT.mkdir(parents=True, exist_ok=True)
        SITES.write_text(json.dumps(known))
    return {n: known[n] for n in names if known.get(n)}


def _km(lat1, lon1, lat2, lon2):
    return 111.2 * np.hypot(lat2 - lat1, (lon2 - lon1) * np.cos(np.radians(lat1)))


def site_for(lat, lon):
    """The nearest MOP site to a break, or None outside California or beyond MAX_KM."""
    if not (BOX[0] <= lat <= BOX[1] and BOX[2] <= lon <= BOX[3]):
        return None
    cat = _catalog()
    coarse = _sites([s for ids in cat.values() for s in ids[::COARSE]])
    best = min(coarse, key=lambda s: _km(lat, lon, *coarse[s][:2]))
    ids = cat[best.rstrip("0123456789")]
    k = ids.index(best)
    near = _sites(ids[max(0, k - COARSE):k + COARSE + 1])
    site = min(near, key=lambda s: _km(lat, lon, *near[s][:2]))
    return site if _km(lat, lon, *near[site][:2]) <= MAX_KM else None


def _path(site):
    return OUT / f"{site}.nc"


def fetch_site(site, lon):
    """Morning spectra at one site, from the WW3 archive's start (2017) to now: hindcast and nowcast."""
    path = _path(site)
    if path.exists():
        return
    hour = int(round(config.LOCAL_HOUR - lon / 15)) % 24
    parts = []
    for kind in ("hindcast", "nowcast"):
        with xr.open_dataset(URL.format(site, kind)) as d:
            t = d.waveTime.values
            keep = np.flatnonzero((t >= np.datetime64("2016-12-31")) & (t.astype("datetime64[h]").astype(int) % 24 == hour))
            if not len(keep):
                continue
            sub = d[["waveHs", "waveTp", "waveDp", "waveEnergyDensity", "waveMeanDirection"]]
            sub = sub.isel(waveTime=slice(keep[0], keep[-1] + 1)).load()
            parts.append(sub.isel(waveTime=keep - keep[0]))
            freq, band = d.waveFrequency.values, d.waveBandwidth.values
            meta = dict(lat=float(d.metaLatitude), lon=float(d.metaLongitude), depth=float(d.metaWaterDepth),
                        shore_normal=float(d.metaShoreNormal))
    ds = xr.concat(parts, "waveTime").assign_coords(waveFrequency=freq)
    ds["waveBandwidth"] = ("waveFrequency", band)
    ds.attrs.update(meta, site=site)
    OUT.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path.with_suffix(".part"))
    path.with_suffix(".part").rename(path)


def fetch():
    """Find and download the MOP site for every break in California."""
    from .points import all_points

    for kind, name, region, lat, lon, facing in all_points():
        site = site_for(lat, lon)
        if site:
            fetch_site(site, lon)
            print(f"  {name}: {site}", flush=True)


def _trains(e, f, df, dirs):
    """Split one spectrum at the troughs between its peaks: [(hs, tp, dir)] biggest first."""
    peaks = [k for k in range(len(e)) if e[k] > 0.05 * e.max()
             and (k == 0 or e[k] >= e[k - 1]) and (k == len(e) - 1 or e[k] > e[k + 1])]
    cuts = []
    for a, b in zip(peaks, peaks[1:]):
        trough = a + int(np.argmin(e[a:b + 1]))
        if e[trough] < 0.75 * min(e[a], e[b]):  # well separated: two trains
            cuts.append(trough)
    out = []
    for lo, hi in zip([0] + cuts, cuts + [len(e)]):
        seg = slice(lo, hi)
        k = lo + int(np.argmax(e[seg]))
        out.append((4 * np.sqrt(np.sum(e[seg] * df[seg])), 1 / f[k], dirs[k]))
    return sorted(out, key=lambda t: -t[0])


def mornings(lat, lon):
    """Wave arrays like points.mornings() from the break's MOP site, or None without one:
    day, hs1/tp1/dir1 and hs2/tp2/dir2 (the two biggest trains), hs/tp/dir (all waves), cell, site."""
    site = site_for(lat, lon)
    if site is None or not _path(site).exists():
        return None
    with xr.open_dataset(_path(site)) as d:
        f, df = d.waveFrequency.values, d.waveBandwidth.values
        e, dirs = d.waveEnergyDensity.values, d.waveMeanDirection.values
        local = d.waveTime.values + np.timedelta64(int(round(lon / 15 * 60)), "m")
        out = dict(day=local.astype("datetime64[D]"), hs=d.waveHs.values.astype(float), tp=d.waveTp.values.astype(float),
                   dir=d.waveDp.values.astype(float), cell=(d.attrs["lat"], d.attrs["lon"]), site=site,
                   depth=d.attrs["depth"])
    n = len(out["day"])
    for k in ("hs1", "tp1", "dir1", "hs2", "tp2", "dir2"):
        out[k] = np.full(n, np.nan)
    for i in range(n):
        if not np.isfinite(e[i]).all():
            continue
        for j, (h, t, dr) in enumerate(_trains(e[i], f, df, dirs[i])[:2], 1):
            out[f"hs{j}"][i], out[f"tp{j}"][i], out[f"dir{j}"][i] = h, t, dr
    return out
