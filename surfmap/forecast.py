"""The coming week's mornings at every marked break, in the shape of the climate data, so the atlas
judges each one by the viewer's own criteria.

Waves: CDIP MOP's nearshore forecast in California; elsewhere Copernicus Marine's MFWAM forecast at
the break's own ~8 km cell (the same model and cell as the ten-year record) when a Copernicus login is
available, else Open-Meteo's marine forecast. Wind: Open-Meteo's 10 m wind and gusts just offshore.
Tide is left to the page, which has each break's tide constituents. Each day is sampled at every
hour in config.HOURS (local solar time), dated by the local day. Water: the latest satellite reading
at each break (water()), GLERL's GLSEA in the Great Lakes and NOAA's OISST elsewhere, for the wetsuit.

Reads the break list from web/data/meta.json and writes web/data/forecast.json:
    {"generated": ISO time, "sources": {...}, "hours": [...], "points": {"name@lat,lon": {"8": [[day,
     utc hour, hs1, tp1, dir1, hs2, tp2, dir2, hs, tp, dir, kt, wind_from, gust_kt], ...], "12": [...], ...}},
     "water": {"name@lat,lon": [°C, day, source]}}

chart() writes the map's swell chart: the whole ocean on a 1 degree grid every CHART_STEP_H hours,
wave height, peak period and direction from NOAA's WAVEWATCH III and 10 m wind from NCEP's GFS,
both from NOAA CoastWatch's ERDDAP (no login). web/data/fcmap.json describes it (times, grid, the
fields and their steps) and web/data/fcmap.b64.txt holds it: gzip, base64, Uint8 [time][field][lat][lon],
255 where there's no value. WAVEWATCH III leaves out the Great Lakes and Hudson Bay, so the enclosed-sea
points of the map (meta.json cells marked enclosed) get Open-Meteo's marine forecast as well, in
fcmap.json "basin": {cells: their indexes, data: gzip base64 Uint8 [time][hs, tp, dir][cell]}.
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
PACE = 6     # seconds between Open-Meteo requests: each location counts, and the free tier allows 600 a minute
COASTWATCH = "https://coastwatch.pfeg.noaa.gov/erddap/griddap"
GLERL = "https://apps.glerl.noaa.gov/erddap/griddap"
CHART_STEP_H = 6
CHART_FIELDS = [("hs", 0.1), ("tp", 0.1), ("dir", 2), ("kt", 0.5), ("wind_from", 2)]  # name, units per byte step
MAP_JSON, MAP_DATA = OUT.parent / "fcmap.json", OUT.parent / "fcmap.b64.txt"


def key(p):
    return f"{p['name']}@{p['lat']},{p['lon']}"


def _now():
    return datetime.now(timezone.utc).replace(tzinfo=None, minute=0, second=0, microsecond=0)


def _rows(times, lon, values, step_h):
    """{hour: {local day: (utc hour, values)}}: the step nearest each config.HOURS solar time, each day."""
    from .points import _local_mornings

    t = np.asarray(times, dtype="datetime64[m]")
    out = {}
    for hour in config.HOURS:
        idx, days = _local_mornings(t, lon, step_h, hour)
        out[hour] = {str(d): (str(np.datetime64(times[i], "h")), [v[i] for v in values]) for i, d in zip(idx, days)}
    return out


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

    from .points import corrected, correction, sea_tp, trains

    la, lo = p["cell"]
    ds = cm.open_dataset(dataset_id="cmems_mod_glo_wav_anfc_0.083deg_PT3H-i",
                         variables=["VHM0", "VTM10", "VMDR", "VHM0_SW1", "VTM01_SW1", "VMDR_SW1", "VHM0_SW2", "VTM01_SW2", "VMDR_SW2"],
                         minimum_longitude=lo - 0.15, maximum_longitude=lo + 0.15, minimum_latitude=la - 0.15, maximum_latitude=la + 0.15,
                         start_datetime=(_now() - timedelta(hours=3)).isoformat(), end_datetime=(_now() + timedelta(days=DAYS)).isoformat()).load()
    wet = np.isfinite(ds.VHM0).any("time").values  # the nearest sea cell to the record's own (which may be on the 0.2° grid)
    d = np.where(wet, np.hypot(*np.meshgrid(ds.longitude - lo, ds.latitude - la)), np.inf)
    j, i = np.unravel_index(np.argmin(d), d.shape)
    c = ds.isel(latitude=j, longitude=i)
    g = lambda v: c[v].values.astype(float)
    t, hs = trains(g("VHM0_SW1"), g("VTM01_SW1"), g("VMDR_SW1"), g("VHM0_SW2"), g("VTM01_SW2"), g("VMDR_SW2")), g("VHM0")
    qm = correction(p["lat"], p["lon"])
    if qm is not None:  # as points.mornings() corrects the record
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.where(hs > 0, corrected(hs, qm) / hs, 1.0)
        t["hs1"], t["hs2"], hs = t["hs1"] * f, t["hs2"] * f, hs * f
    cols = [t["hs1"], t["tp1"], t["dir1"], t["hs2"], t["tp2"], t["dir2"], hs, sea_tp(g("VTM10"), p["measure"] == "waves"), g("VMDR")]
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
        time.sleep(PACE)
    return [b.get("hourly") for b in out]


def _marine_fallback(pts):
    """Open-Meteo's marine forecast for points without CDIP or Copernicus waves. Its swell peak
    periods are often missing, so both trains go from their mean periods through the same
    conversion as the record (points.trains)."""
    from .points import corrected, correction, sea_tp, trains

    names = ["swell_wave_height", "swell_wave_period", "swell_wave_direction", "secondary_swell_wave_height",
             "secondary_swell_wave_period", "secondary_swell_wave_direction", "wave_height", "wave_period", "wave_direction"]
    out = {}
    for p, h in zip(pts, _open_meteo(MARINE, [(p["lat"], p["lon"]) for p in pts], names)):
        if not h:
            continue
        cols = [np.array([np.nan if v is None else v for v in h[n]], float) for n in names]
        t, hs = trains(*cols[:6]), cols[6]
        qm = correction(p["lat"], p["lon"])
        if qm is not None:  # as points.mornings() corrects the record
            with np.errstate(invalid="ignore", divide="ignore"):
                f = np.where(hs > 0, corrected(hs, qm) / hs, 1.0)
            t["hs1"], t["hs2"], hs = t["hs1"] * f, t["hs2"] * f, hs * f
        cols = [t["hs1"], t["tp1"], t["dir1"], t["hs2"], t["tp2"], t["dir2"], hs, sea_tp(cols[7], p["measure"] == "waves"), cols[8]]
        out[key(p)] = _rows(np.array(h["time"], dtype="datetime64[m]"), p["lon"], cols, 1)
    return out


def _wind(pts):
    from .points import WIND_DEG, _out

    spots = [tuple(float(x) for x in _out(p["lat"], p["lon"], p["facing"], WIND_DEG)) for p in pts]
    out = {}
    names = ["wind_speed_10m", "wind_direction_10m", "wind_gusts_10m"]
    for p, h in zip(pts, _open_meteo(WEATHER, spots, names, dict(wind_speed_unit="kn"))):
        if h:
            cols = [np.array([np.nan if v is None else v for v in h[n]], float) for n in names]
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
    try:
        temps = water(pts)
    except Exception as e:  # the page falls back on the monthly means
        print(f"  water temperature left out: {type(e).__name__}: {str(e)[:120]}", flush=True)
        temps = {}

    out = {}
    for p in pts:
        k = key(p)
        per = {}
        for h in config.HOURS:
            w, v = waves.get(k, {}).get(h, {}), wind.get(k, {}).get(h, {})
            rows = [[day, utc, *[None if not np.isfinite(x) else round(float(x), 2) for x in vals],
                     *[None if not np.isfinite(x) else round(float(x), 1) for x in v.get(day, (None, [np.nan] * 3))[1]]]
                    for day, (utc, vals) in sorted(w.items())]
            if rows:
                per[str(h)] = rows[:DAYS]
        if per:
            out[k] = per
    OUT.write_text(json.dumps(dict(generated=_now().isoformat() + "Z", sources=sources, hours=list(config.HOURS), points=out, water=temps),
                              separators=(",", ":")))
    print(f"forecast: {len(out)} of {len(pts)} breaks; waves from {sources}; wrote {OUT}", flush=True)


def _erddap(dataset, variables, sel, server=COASTWATCH):
    """One griddap subset (CoastWatch by default) as a loaded xarray Dataset."""
    import tempfile

    import xarray as xr

    from .ww3 import _quote

    url = f"{server}/{dataset}.nc?" + ",".join(_quote(v + sel) for v in variables)
    for attempt in range(4):
        r = requests.get(url, timeout=900)
        if r.ok:
            break
        time.sleep(60 * (attempt + 1))
    r.raise_for_status()
    with tempfile.NamedTemporaryFile(suffix=".nc") as f:
        f.write(r.content)
        f.flush()
        with xr.open_dataset(f.name) as d:
            return d.load()


def chart():
    """The coming week's swell chart for the map (see the module notes)."""
    import base64
    import gzip

    from .points import KT

    now = _now()
    t0 = now.replace(hour=now.hour - now.hour % CHART_STEP_H)
    stamp = f"{t0:%Y-%m-%dT%H}:00:00Z"
    waves = _erddap("NWW3_Global_Best", ["Thgt", "Tper", "Tdir"],
                    f"[({stamp}):{CHART_STEP_H}:last][(0.0)][(-77.5):2:(77.5)][(0.0):2:(359.5)]").squeeze("depth", drop=True)
    wind = _erddap("NCEP_Global_Best", ["ugrd10m", "vgrd10m"], f"[({stamp}):{CHART_STEP_H // 3}:last][(77.5):2:(-77.5)][(0.0):2:(359.5)]")
    wind = wind.sortby("latitude").sel(time=waves.time, method="nearest", tolerance=np.timedelta64(90, "m"))
    u, v = wind.ugrd10m.values, wind.vgrd10m.values
    fields = {"hs": waves.Thgt.values, "tp": waves.Tper.values, "dir": waves.Tdir.values % 360,
              "kt": np.hypot(u, v) * KT, "wind_from": np.degrees(np.arctan2(-u, -v)) % 360}
    with np.errstate(invalid="ignore"):
        cube = np.stack([np.where(np.isfinite(fields[n]), np.clip(np.round(fields[n] / step), 0, 254), 255)
                         for n, step in CHART_FIELDS], axis=1).astype(np.uint8)  # [time][field][lat][lon]
    lat, lon = waves.latitude.values, waves.longitude.values
    MAP_DATA.write_text(base64.b64encode(gzip.compress(cube.tobytes(), 9)).decode())
    try:
        basin = _chart_basins(waves.time.values)
    except Exception as e:  # rate limits, outages: the chart still has the open ocean
        print(f"  enclosed seas left out: {type(e).__name__}: {str(e)[:120]}", flush=True)
        basin = None
    MAP_JSON.write_text(json.dumps(dict(
        generated=now.isoformat() + "Z", sources={"waves": "NOAA WAVEWATCH III (NWW3_Global_Best)", "wind": "NCEP GFS (NCEP_Global_Best)",
                                                  "enclosed seas": "Open-Meteo marine forecast"},
        times=[str(np.datetime64(t, "m")) + "Z" for t in waves.time.values], lat0=float(lat[0]), lon0=float(lon[0]),
        step=float(lat[1] - lat[0]), nlat=len(lat), nlon=len(lon), fields=CHART_FIELDS, basin=basin), separators=(",", ":")))
    print(f"swell chart: {cube.shape[0]} times, {len(lat)} x {len(lon)} grid; wrote {MAP_DATA} "
          f"({MAP_DATA.stat().st_size / 1e6:.1f} MB)", flush=True)


def _chart_basins(times):
    """Open-Meteo's waves at the map's enclosed-sea points at the chart's times (see chart())."""
    import base64
    import gzip

    from .points import corrected, correction, sea_tp

    cells = json.loads(META.read_text())["cells"]
    idx = [k for k, e in enumerate(cells["enclosed"]) if e]
    at = [(cells["lat"][k], cells["lon"][k]) for k in idx]
    want = np.asarray(times, dtype="datetime64[m]")
    cube = np.full((len(want), 3, len(idx)), 255, np.uint8)
    for n, h in enumerate(_open_meteo(MARINE, at, ["wave_height", "wave_period", "wave_direction"])):
        if not h:
            continue
        t = np.array(h["time"], dtype="datetime64[m]")
        pos = np.minimum(np.searchsorted(t, want), len(t) - 1)
        hit = t[pos] == want
        hs, tm, d = (np.array([np.nan if v is None else v for v in h[k]], float)[pos] for k in ("wave_height", "wave_period", "wave_direction"))
        qm = correction(*at[n])
        if qm is not None:
            hs = corrected(hs, qm)
        for f, (v, step) in enumerate([(hs, 0.1), (sea_tp(tm, True), 0.1), (d % 360, 2)]):
            ok = hit & np.isfinite(v)
            cube[ok, f, n] = np.clip(np.round(v[ok] / step), 0, 254)
    print(f"  enclosed seas: {len(idx)} points", flush=True)
    return dict(cells=idx, data=base64.b64encode(gzip.compress(cube.tobytes(), 9)).decode())


def water(pts):
    """{key: [°C, day, source]}: the latest satellite water temperature at each break, at the nearest
    wet cell within half a degree: GLERL's GLSEA (~1.4 km) in the Great Lakes, NOAA's OISST (0.25°) elsewhere."""
    from .points import lake_of

    grids = {"OISST": _erddap("ncdcOisst21NrtAgg_LonPM180", ["sst"], "[last][0][(-89.875):1:(89.875)][(-179.875):1:(179.875)]").squeeze("zlev", drop=True)}
    if any(lake_of(p["lat"], p["lon"]) for p in pts):
        grids["GLSEA"] = _erddap("GLSEA_ACSPO_GCS", ["sst"], "[last][(38.87):1:(50.61)][(-92.42):1:(-75.88)]", server=GLERL)
    out = {}
    for p in pts:
        src = "GLSEA" if "GLSEA" in grids and lake_of(p["lat"], p["lon"]) else "OISST"
        g = grids[src]
        la, lo, v = g.latitude.values, g.longitude.values, g.sst.values[0]
        i, j = np.flatnonzero(np.abs(la - p["lat"]) <= 0.5), np.flatnonzero(np.abs((lo - p["lon"] + 180) % 360 - 180) <= 0.5)
        box = v[np.ix_(i, j)]
        if not np.isfinite(box).any():
            continue
        d = np.where(np.isfinite(box), np.hypot(*np.meshgrid(lo[j] - p["lon"], la[i] - p["lat"])[::-1]), np.inf)
        a, b = np.unravel_index(np.argmin(d), d.shape)
        out[key(p)] = [round(float(box[a, b]), 1), str(g.time.values[0])[:10], src]
    print(f"  water temperature at {len(out)} breaks", flush=True)
    return out
