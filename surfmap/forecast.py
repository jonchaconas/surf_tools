"""The coming week's mornings at every marked break, in the shape of the climate data, so the atlas
judges each one by the viewer's own criteria.

Waves: CDIP MOP's nearshore forecast in California; elsewhere Copernicus Marine's MFWAM forecast at
the break's own ~8 km cell (the same model and cell as the ten-year record) when a Copernicus login is
available, else Open-Meteo's marine forecast. Wind: Open-Meteo's 10 m wind forecast just offshore.
Tide is left to the page, which has each break's tide constituents. Each morning is the step
nearest config.LOCAL_HOUR solar time, dated by the local day.

Reads the break list from web/data/meta.json and writes web/data/forecast.json:
    {"generated": ISO time, "sources": {...}, "points": {"name@lat,lon": [[day, utc hour,
     hs1, tp1, dir1, hs2, tp2, dir2, hs, tp, dir, kt, wind_from], ...]}}
"""

import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import numpy as np
import requests

from . import config

META = config.ROOT / "web" / "data" / "meta.json"
OUT = config.ROOT / "web" / "data" / "forecast.json"
DAYS = 8
MARINE = "https://marine-api.open-meteo.com/v1/marine"
WEATHER = "https://api.open-meteo.com/v1/forecast"
BATCH = 50  # locations per Open-Meteo request


def key(p):
    return f"{p['name']}@{p['lat']},{p['lon']}"


def _now():
    return datetime.now(timezone.utc).replace(tzinfo=None, minute=0, second=0, microsecond=0)


def _rows(times, lon, values, step_h):
    """{local day: (utc hour, values at the morning step)} for the step nearest LOCAL_HOUR each day."""
    from .points import _local_mornings

    idx, days = _local_mornings(np.asarray(times, dtype="datetime64[m]"), lon, step_h)
    return {str(d): (str(np.datetime64(times[i], "h")), [v[i] for v in values]) for i, d in zip(idx, days)}


def _cdip(p):
    """California: CDIP's nearshore forecast, the two biggest trains from each morning's spectrum."""
    import xarray as xr

    from .cdip import URL, _trains

    with xr.open_dataset(URL.format(p["cdip"], "forecast")) as d:
        t = d.waveTime.values
        keep = t >= np.datetime64(_now() - timedelta(hours=3))
        e, dirs = d.waveEnergyDensity.values[keep], d.waveMeanDirection.values[keep]
        f, df = d.waveFrequency.values, d.waveBandwidth.values
        hs, tp, dp = d.waveHs.values[keep], d.waveTp.values[keep], d.waveDp.values[keep]
        t = t[keep]
    trains = [_trains(e[i], f, df, dirs[i])[:2] + [(np.nan, np.nan, np.nan)] * 2 for i in range(len(t))]
    cols = [[tr[j][k] for tr in trains] for j in range(2) for k in range(3)] + [hs, tp, dp]
    return _rows(t, p["lon"], cols, 3)


def _copernicus(p):
    """MFWAM's forecast at the break's own cell (meta 'cell'), both swell partitions."""
    import copernicusmarine as cm

    from .inland import TP_RATIO

    la, lo = p["cell"]
    ds = cm.open_dataset(dataset_id="cmems_mod_glo_wav_anfc_0.083deg_PT3H-i",
                         variables=["VHM0", "VTM10", "VMDR", "VHM0_SW1", "VTM01_SW1", "VMDR_SW1", "VHM0_SW2", "VTM01_SW2", "VMDR_SW2"],
                         minimum_longitude=lo - 0.02, maximum_longitude=lo + 0.02, minimum_latitude=la - 0.02, maximum_latitude=la + 0.02,
                         start_datetime=(_now() - timedelta(hours=3)).isoformat(), end_datetime=(_now() + timedelta(days=DAYS)).isoformat())
    c = ds.sel(latitude=la, longitude=lo, method="nearest").load()
    g = lambda v: c[v].values.astype(float)
    cols = [g("VHM0_SW1"), g("VTM01_SW1") * TP_RATIO, g("VMDR_SW1"), g("VHM0_SW2"), g("VTM01_SW2") * TP_RATIO, g("VMDR_SW2"),
            g("VHM0"), g("VTM10") * TP_RATIO, g("VMDR")]
    return _rows(c.time.values, p["lon"], cols, 3)


def _open_meteo(url, pts, hourly, extra=None):
    """One hourly series per point, in batches; None for a point the API has nothing for."""
    out = []
    for i in range(0, len(pts), BATCH):
        chunk = pts[i:i + BATCH]
        params = dict(latitude=",".join(f"{a:.3f}" for a, _ in chunk), longitude=",".join(f"{b:.3f}" for _, b in chunk),
                      hourly=",".join(hourly), timezone="GMT", forecast_days=DAYS, **(extra or {}))
        for attempt in range(4):
            r = requests.get(url, params=params, timeout=120)
            if r.status_code == 429:
                time.sleep(60 * (attempt + 1))
                continue
            r.raise_for_status()
            break
        body = r.json()
        out += body if isinstance(body, list) else [body]
    return [b.get("hourly") for b in out]


def _marine_fallback(pts):
    """Open-Meteo's marine forecast for points without CDIP or Copernicus waves."""
    names = ["swell_wave_height", "swell_wave_peak_period", "swell_wave_direction", "secondary_swell_wave_height",
             "secondary_swell_wave_period", "secondary_swell_wave_direction", "wave_height", "wave_period", "wave_direction"]
    out = {}
    for p, h in zip(pts, _open_meteo(MARINE, [(p["lat"], p["lon"]) for p in pts], names)):
        if not h:
            continue
        cols = [np.array([np.nan if v is None else v for v in h[n]], float) for n in names]
        cols[7] = cols[7] * 1.29  # mean -> peak period, as elsewhere in the repo
        out[key(p)] = _rows(np.array(h["time"], dtype="datetime64[m]"), p["lon"], cols, 1)
    return out


def _wind(pts):
    from .points import WIND_DEG, _out

    spots = [tuple(float(x) for x in _out(p["lat"], p["lon"], p["facing"], WIND_DEG)) for p in pts]
    out = {}
    for p, h in zip(pts, _open_meteo(WEATHER, spots, ["wind_speed_10m", "wind_direction_10m"], dict(wind_speed_unit="kn"))):
        if h:
            cols = [np.array([np.nan if v is None else v for v in h[n]], float) for n in ("wind_speed_10m", "wind_direction_10m")]
            out[key(p)] = _rows(np.array(h["time"], dtype="datetime64[m]"), p["lon"], cols, 1)
    return out


def build():
    pts = json.loads(META.read_text())["points"]
    waves, sources = {}, {"cdip": 0, "copernicus": 0, "open-meteo": 0}

    def one(p):
        try:
            if p.get("cdip"):
                return key(p), "cdip", _cdip(p)
            return key(p), "copernicus", _copernicus(p)
        except Exception as e:  # no login, outage: Open-Meteo below
            print(f"  {p['name']}: {type(e).__name__}: {str(e)[:120]}", flush=True)
            return key(p), None, None

    with ThreadPoolExecutor(6) as pool:
        for k, src, rows in pool.map(one, pts):
            if rows:
                waves[k] = rows
                sources[src] += 1
    rest = [p for p in pts if key(p) not in waves]
    if rest:
        waves.update(_marine_fallback(rest))
        sources["open-meteo"] = sum(key(p) in waves for p in rest)
    wind = _wind(pts)

    out = {}
    for p in pts:
        k = key(p)
        w, v = waves.get(k, {}), wind.get(k, {})
        rows = [[day, hour, *[None if not np.isfinite(x) else round(float(x), 2) for x in vals],
                 *[None if not np.isfinite(x) else round(float(x), 1) for x in v.get(day, (None, [np.nan, np.nan]))[1]]]
                for day, (hour, vals) in sorted(w.items())]
        if rows:
            out[k] = rows[:DAYS]
    OUT.write_text(json.dumps(dict(generated=_now().isoformat() + "Z", sources=sources, points=out), separators=(",", ":")))
    print(f"forecast: {len(out)} of {len(pts)} breaks; waves from {sources}; wrote {OUT}", flush=True)
