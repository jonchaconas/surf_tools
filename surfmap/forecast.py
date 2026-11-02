"""The coming week's mornings at every marked break, in the shape of the climate data, so the atlas
judges each one by the viewer's own criteria.

Waves: NOAA's Great Lakes wave model (GLWU, 2.5 km) on the Great Lakes; CDIP MOP's nearshore forecast in
California; elsewhere Copernicus Marine's MFWAM forecast at
the break's own ~8 km cell (the same model and cell as the ten-year record) when a Copernicus login is
available, else Open-Meteo's marine forecast. Weather: Open-Meteo's 10 m wind and gusts, air temperature
and rain just offshore (GFS from the swell chart's download when Open-Meteo is out of quota).
Tide is left to the page, which has each break's tide constituents. Each day is sampled at every
hour in config.HOURS (local solar time), dated by the local day. Water: the latest satellite reading
at each break (water()), GLERL's GLSEA in the Great Lakes and NOAA's OISST elsewhere, for the wetsuit.

Reads the break list from web/data/meta.json and writes web/data/forecast.json:
    {"generated": ISO time, "sources": {...}, "hours": [...], "points": {"name@lat,lon": {"8": [[day,
     utc hour, hs1, tp1, dir1, hs2, tp2, dir2, hs, tp, dir, kt, wind_from, gust_kt, air_c, rain_mm,
     rain_chance_pct], ...], "12": [...], ...}},
     "water": {"name@lat,lon": [°C, day, source]}}

chart() writes the map's swell chart: the whole ocean on a 1 degree grid every CHART_STEP_H hours,
wave height, peak period and direction from NOAA's GFS-Wave (NOMADS; WAVEWATCH III via CoastWatch if
that fails) and 10 m wind from NCEP's GFS (CoastWatch), no login. Enclosed seas also come at their
model's own grid (CHART_REGIONS: GLWU 0.1° for the Great Lakes, GFS-Wave 0.25° for the rest), in
fcmap.json "regions": [{name, lat0, lon0, step, nlat, nlon, data: gzip base64 Uint8 [time][hs, tp, dir][lat][lon]}]. web/data/fcmap.json describes it (times, grid, the
fields and their steps) and web/data/fcmap.b64.txt holds it: gzip, base64, Uint8 [time][field][lat][lon],
255 where there's no value. WAVEWATCH III leaves out the Great Lakes and Hudson Bay, so the enclosed-sea
points of the map (meta.json cells marked enclosed) get Open-Meteo's marine forecast as well, in
fcmap.json "basin": {cells: their indexes, data: gzip base64 Uint8 [time][hs, tp, dir][cell]}.
"""

import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from functools import lru_cache

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
GLWU = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/glwu/prod"  # NOAA's Great Lakes wave model, GRIB2
GFSWAVE = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"  # NOAA's global wave model, GRIB2
CHART_FRAMES = 48  # six days of 3-hour steps
# Seas the 1° chart is too coarse for, sent at their model's own grid: (lat0, lat1, lon0, lon1, source).
# The Black Sea and the Caspian are in neither NOAA model; they keep the enclosed-sea points (Open-Meteo).
CHART_REGIONS = {
    "Great Lakes": (41.0, 49.5, -92.6, -75.8, "glwu"),
    "Mediterranean": (30.0, 46.0, -6.0, 37.0, "gfswave"),
    "Baltic": (53.5, 66.0, 9.0, 31.0, "gfswave"),
    "Red Sea": (12.0, 30.0, 32.0, 44.0, "gfswave"),
    "Persian Gulf": (23.0, 31.0, 47.0, 57.0, "gfswave"),
    "Hudson Bay": (50.5, 65.0, -96.0, -76.0, "gfswave"),
}
REGION_STEP = {"glwu": 0.1, "gfswave": 0.25}
GLWU_FIELDS = {"HTSGW": "hs", "PERPW": "tp", "DIRPW": "dir", "WIND": "kt", "WDIR": "wind_from"}
CHART_STEP_H = 3  # GFS's own step; WAVEWATCH III is hourly
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
            if r.status_code == 429 and "Daily" not in r.text:  # the minute's quota clears; the day's doesn't
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


WEATHER_NAMES = ["wind_speed_10m", "wind_direction_10m", "wind_gusts_10m", "temperature_2m", "precipitation", "precipitation_probability"]
_GFS = {}  # air temperature and rain grids from chart()'s GFS download, for _from_chart()


def _weather(pts):
    """{key: rows of wind (kt), its direction, gusts, air °C, rain mm in the hour and its chance %}."""
    from .points import WIND_DEG, _out

    spots = [tuple(float(x) for x in _out(p["lat"], p["lon"], p["facing"], WIND_DEG)) for p in pts]
    out = {}
    names = WEATHER_NAMES
    for p, h in zip(pts, _open_meteo(WEATHER, spots, names, dict(wind_speed_unit="kn"))):
        if h:
            cols = [np.array([np.nan if v is None else v for v in h[n]], float) for n in names]
            out[key(p)] = _rows(np.array(h["time"], dtype="datetime64[m]"), p["lon"], cols, 1)
    return out


def build():
    pts = json.loads(META.read_text())["points"]
    waves, sources = {}, {"glwu": 0, "cdip": 0, "copernicus": 0, "open-meteo": 0, "swell chart": 0}
    waves.update(_glwu_rows(pts))  # the Great Lakes' own model first
    sources["glwu"] = len(waves)

    def one(p):
        try:
            if p.get("cdip"):
                return key(p), "cdip", _cdip(p)
            return key(p), "copernicus", _copernicus(p)
        except Exception as e:  # no login, outage: Open-Meteo below
            print(f"  {p['name']}: {type(e).__name__}: {str(e)[:120]}", flush=True)
            return key(p), None, None

    with ThreadPoolExecutor(6) as pool:
        for k, src, rows in pool.map(one, [p for p in pts if key(p) not in waves]):
            if rows:
                waves[k] = rows
                sources[src] += 1
    rest = [p for p in pts if key(p) not in waves]
    if rest:
        try:
            waves.update(_marine_fallback(rest))
            sources["open-meteo"] = sum(key(p) in waves for p in rest)
        except Exception as e:  # out of quota: the swell chart below
            print(f"  Open-Meteo waves failed: {type(e).__name__}: {str(e)[:120]}", flush=True)
    try:
        wind = _weather(pts)
    except Exception as e:
        print(f"  Open-Meteo wind failed: {type(e).__name__}: {str(e)[:120]}", flush=True)
        wind = {}
    if MAP_JSON.exists() and (len(waves) < len(pts) or len(wind) < len(pts)):
        chart_waves, chart_wind = _from_chart(pts)
        for p in pts:
            if key(p) not in waves and key(p) in chart_waves:
                waves[key(p)] = chart_waves[key(p)]
                sources["swell chart"] += 1
            wind.setdefault(key(p), chart_wind.get(key(p), {}))
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
                     *[None if not np.isfinite(x) else round(float(x), 1) for x in v.get(day, (None, [np.nan] * len(WEATHER_NAMES)))[1]]]
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
    times = np.array([np.datetime64(t0) + np.timedelta64(CHART_STEP_H * k, "h") for k in range(CHART_FRAMES)], dtype="datetime64[m]")
    lat, lon = np.arange(-77.5, 77.6, 1.0), np.arange(0.0, 360.0, 1.0)
    waves, regions = _gfswave_chart(times, lat, lon)
    if waves is None:  # NOMADS down: CoastWatch's WAVEWATCH III, which leaves out most enclosed seas
        stamp = f"{t0:%Y-%m-%dT%H}:00:00Z"
        w = _erddap("NWW3_Global_Best", ["Thgt", "Tper", "Tdir"],
                    f"[({stamp}):{CHART_STEP_H}:last][(0.0)][(-77.5):2:(77.5)][(0.0):2:(359.5)]").squeeze("depth", drop=True)
        times = w.time.values.astype("datetime64[m]")
        waves = {"hs": w.Thgt.values, "tp": w.Tper.values, "dir": w.Tdir.values % 360}
    stamp = f"{t0:%Y-%m-%dT%H}:00:00Z"
    wind = _erddap("NCEP_Global_Best", ["ugrd10m", "vgrd10m", "tmp2m", "pratesfc"], f"[({stamp}):{CHART_STEP_H // 3}:last][(77.5):2:(-77.5)][(0.0):2:(359.5)]")
    wind = wind.sortby("latitude").sel(time=times, method="nearest", tolerance=np.timedelta64(90, "m"))
    u, v = wind.ugrd10m.values, wind.vgrd10m.values
    _GFS.update(air=wind.tmp2m.values - 273.15, rain=wind.pratesfc.values * 3600)  # K -> °C, kg m-2 s-1 -> mm an hour
    fields = {**waves, "kt": np.hypot(u, v) * KT, "wind_from": np.degrees(np.arctan2(-u, -v)) % 360}
    cube = np.stack([_bytes(fields[n], step) for n, step in CHART_FIELDS], axis=1)  # [time][field][lat][lon]
    MAP_DATA.write_text(base64.b64encode(gzip.compress(cube.tobytes(), 9)).decode())
    regions += [r for r in [_glwu_region(times)] if r]
    try:
        basin = _chart_basins(times)
    except Exception as e:  # rate limits, outages: the chart still has the open ocean
        print(f"  enclosed seas left out: {type(e).__name__}: {str(e)[:120]}", flush=True)
        basin = None
    MAP_JSON.write_text(json.dumps(dict(
        generated=now.isoformat() + "Z", sources={"waves": "NOAA GFS-Wave (NOMADS)", "wind": "NCEP GFS (NCEP_Global_Best)",
                                                  "Great Lakes": "NOAA GLWU", "Black Sea, Caspian": "Open-Meteo marine forecast"},
        times=[str(t) + "Z" for t in times], lat0=float(lat[0]), lon0=float(lon[0]), step=1.0, nlat=len(lat), nlon=len(lon),
        fields=CHART_FIELDS, basin=basin, regions=regions), separators=(",", ":")))
    print(f"swell chart: {cube.shape[0]} times, {len(lat)} x {len(lon)} grid, regions {[r['name'] for r in regions]}; "
          f"wrote {MAP_DATA} ({MAP_DATA.stat().st_size / 1e6:.1f} MB)", flush=True)


def _bytes(v, step):
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(v), np.clip(np.round(v / step), 0, 254), 255).astype(np.uint8)


def _region_dict(name, lat0, lon0, step, cube):
    import base64
    import gzip

    return dict(name=name, lat0=round(float(lat0), 4), lon0=round(float(lon0), 4), step=step, nlat=cube.shape[2], nlon=cube.shape[3],
                data=base64.b64encode(gzip.compress(cube.tobytes(), 9)).decode())


def _gfswave_chart(times, lat, lon):
    """GFS-Wave's latest run at the chart's times: ({hs, tp, dir} on the 1° chart grid [time, lat, lon],
    the CHART_REGIONS it covers at 0.25°), or (None, []) if NOMADS has no run. Only the three fields'
    bytes are fetched from each forecast hour's file, by its index; one file at a time, so memory stays small."""
    import re
    import tempfile

    import eccodes

    s = requests.Session()
    s.headers["User-Agent"] = "surf_tools (https://github.com/jonchaconas/surf_tools)"
    try:
        run = None
        for day in sorted(set(re.findall(r"gfs\.(\d{8})/", s.get(f"{GFSWAVE}/", timeout=60).text)), reverse=True):
            for cycle in sorted(set(re.findall(r'href="(\d{2})/"', s.get(f"{GFSWAVE}/gfs.{day}/", timeout=60).text)), reverse=True):
                cyc_t = np.datetime64(f"{day[:4]}-{day[4:6]}-{day[6:]}T{cycle}:00")
                last = int((times[-1] - cyc_t) / np.timedelta64(1, "h"))
                if s.head(f"{GFSWAVE}/gfs.{day}/{cycle}/wave/gridded/gfswave.t{cycle}z.global.0p25.f{last:03d}.grib2.idx", timeout=60).ok:
                    run = (day, cycle, cyc_t)
                    break
            if run:
                break
        if run is None or (times[0] - run[2]) < np.timedelta64(0, "h"):
            return None, []
        day, cycle, cyc_t = run
        out = {n: np.full((len(times), len(lat), len(lon)), np.nan, np.float32) for n in ("hs", "tp", "dir")}
        boxes = {name: b for name, b in CHART_REGIONS.items() if b[4] == "gfswave"}
        cubes, axes = {}, None
        for k, t in enumerate(times):
            fh = int((t - cyc_t) / np.timedelta64(1, "h"))
            url = f"{GFSWAVE}/gfs.{day}/{cycle}/wave/gridded/gfswave.t{cycle}z.global.0p25.f{fh:03d}.grib2"
            idx = [line.split(":") for line in s.get(url + ".idx", timeout=60).text.splitlines()]
            want = [i for i, e in enumerate(idx) if e[3] in ("HTSGW", "PERPW", "DIRPW")]
            end = int(idx[want[-1] + 1][1]) - 1 if want[-1] + 1 < len(idx) else ""
            r = s.get(url, headers={"Range": f"bytes={idx[want[0]][1]}-{end}"}, timeout=120)
            r.raise_for_status()
            with tempfile.TemporaryFile() as f:
                f.write(r.content)
                f.seek(0)
                for _ in range(want[-1] - want[0] + 1):
                    h = eccodes.codes_grib_new_from_file(f)
                    try:
                        name = {"swh": "hs", "perpw": "tp", "dirpw": "dir"}.get(eccodes.codes_get(h, "shortName"))
                        if axes is None:
                            nj, ni = eccodes.codes_get(h, "Nj"), eccodes.codes_get(h, "Ni")
                            g_lat = eccodes.codes_get_array(h, "latitudes").reshape(nj, ni)[:, 0]
                            g_lon = eccodes.codes_get_array(h, "longitudes").reshape(nj, ni)[0] % 360
                            axes = (nj, ni, g_lat, g_lon)
                        if name is None:
                            continue
                        v = eccodes.codes_get_values(h)
                        v = np.where(v == eccodes.codes_get(h, "missingValue"), np.nan, v).reshape(axes[0], axes[1])
                    finally:
                        eccodes.codes_release(h)
                    nj, ni, g_lat, g_lon = axes
                    rows = np.searchsorted(-g_lat, -lat)  # the 0.25° grid runs north to south
                    out[name][k] = v[np.ix_(rows, np.searchsorted(g_lon, lon))]
                    for rname, (la0, la1, lo0, lo1, _) in boxes.items():
                        ri = np.flatnonzero((g_lat >= la0) & (g_lat <= la1))[::-1]  # south to north
                        ci = np.flatnonzero(((g_lon - lo0 % 360) % 360) <= (lo1 - lo0))
                        ci = ci[np.argsort((g_lon[ci] - lo0 % 360) % 360)]
                        if rname not in cubes:
                            cubes[rname] = (g_lat[ri[0]], ((g_lon[ci[0]] + 180) % 360) - 180, np.full((len(times), 3, len(ri), len(ci)), 255, np.uint8))
                        cubes[rname][2][k, ["hs", "tp", "dir"].index(name)] = _bytes(v[np.ix_(ri, ci)], {"hs": 0.1, "tp": 0.1, "dir": 2}[name])
            time.sleep(0.3)  # NOMADS allows about 120 requests a minute
        out["dir"] = out["dir"] % 360
        regions = [_region_dict(name, la, lo, 0.25, cube) for name, (la, lo, cube) in cubes.items()]
        print(f"  GFS-Wave run {day} {cycle}z: {len(times)} steps", flush=True)
        return out, regions
    except Exception as e:
        print(f"  GFS-Wave failed: {type(e).__name__}: {str(e)[:150]}", flush=True)
        return None, []


def _glwu_region(times):
    """GLWU on a regular 0.1° grid over the Great Lakes at the chart's times, nearest of its water points within 0.06°."""
    from scipy.spatial import cKDTree

    g = _glwu()
    if g is None:
        return None
    la0, la1, lo0, lo1, _ = CHART_REGIONS["Great Lakes"]
    step = REGION_STEP["glwu"]
    lats, lons = np.arange(la0, la1 + 1e-6, step), np.arange(lo0, lo1 + 1e-6, step)
    gl, go = np.meshgrid(lats, lons, indexing="ij")
    c = np.cos(np.radians(45))
    d, near = cKDTree(np.c_[g["lat"], g["lon"] * c]).query(np.c_[gl.ravel(), go.ravel() * c])
    wet = d <= 0.06
    want = np.asarray(times, dtype="datetime64[m]")
    pos = np.minimum(np.searchsorted(g["times"], want), len(g["times"]) - 1)
    hit = g["times"][pos] == want
    cube = np.full((len(want), 3, len(lats), len(lons)), 255, np.uint8)
    for f, (name, step_v) in enumerate([("hs", 0.1), ("tp", 0.1), ("dir", 2)]):
        vals = g[name][pos][:, near]  # [time, grid cell]
        vals = np.where(wet[None, :] & hit[:, None], vals, np.nan)
        cube[:, f] = _bytes(vals % 360 if name == "dir" else vals, step_v).reshape(len(want), len(lats), len(lons))
    return _region_dict("Great Lakes", lats[0], lons[0], step, cube)


def _chart_basins(times):
    """Waves at the map's enclosed-sea points at the chart's times (see chart()): GLWU on the Great
    Lakes, Open-Meteo's marine forecast elsewhere. None if neither has anything."""
    import base64
    import gzip

    from .points import corrected, correction, lake_of, sea_tp

    cells = json.loads(META.read_text())["cells"]
    idx = [k for k, e in enumerate(cells["enclosed"]) if e]
    at = [(cells["lat"][k], cells["lon"][k]) for k in idx]
    want = np.asarray(times, dtype="datetime64[m]")
    cube = np.full((len(want), 3, len(idx)), 255, np.uint8)

    def put(n, hit, hs, tp, d):
        for f, (v, step) in enumerate([(hs, 0.1), (tp, 0.1), (d % 360, 2)]):
            ok = hit & np.isfinite(v)
            cube[ok, f, n] = np.clip(np.round(v[ok] / step), 0, 254)

    try:
        for n, h in enumerate(_open_meteo(MARINE, at, ["wave_height", "wave_period", "wave_direction"])):
            if not h:
                continue
            t = np.array(h["time"], dtype="datetime64[m]")
            pos = np.minimum(np.searchsorted(t, want), len(t) - 1)
            hs, tm, d = (np.array([np.nan if v is None else v for v in h[k]], float)[pos] for k in ("wave_height", "wave_period", "wave_direction"))
            qm = correction(*at[n])
            put(n, t[pos] == want, corrected(hs, qm) if qm is not None else hs, sea_tp(tm, True), d)
    except Exception as e:  # out of quota: the Great Lakes still come from GLWU below
        print(f"  Open-Meteo enclosed seas failed: {type(e).__name__}: {str(e)[:100]}", flush=True)
    g = _glwu()
    lakes = [n for n, (la, lo) in enumerate(at) if lake_of(la, lo)]
    if g is not None and lakes:
        pos = np.minimum(np.searchsorted(g["times"], want), len(g["times"]) - 1)
        hit = g["times"][pos] == want
        for n in lakes:
            k = _glwu_point(g, *at[n])
            if k is not None:
                put(n, hit, g["hs"][pos, k], g["tp"][pos, k], g["dir"][pos, k])
    filled = (cube != 255).any()
    print(f"  enclosed seas: {len(idx)} points, {'some' if filled else 'none'} filled", flush=True)
    return dict(cells=idx, data=base64.b64encode(gzip.compress(cube.tobytes(), 9)).decode()) if filled else None


@lru_cache(maxsize=1)
def _glwu():
    """NOAA's Great Lakes wave model (GLWU, 2.5 km), its latest run at the chart's steps (CHART_STEP_H
    hours, on the UTC clock): dict(times, lat, lon of its water points, then hs, tp, dir, kt and
    wind_from as [time, point]), or None. Only those fields' bytes are fetched, using the run's index."""
    import re
    import tempfile

    import eccodes

    from .points import KT

    s = requests.Session()
    s.headers["User-Agent"] = "surf_tools (https://github.com/jonchaconas/surf_tools)"
    try:
        days = sorted(set(re.findall(r"glwu\.(\d{8})/", s.get(f"{GLWU}/", timeout=60).text)), reverse=True)
        for day in days:
            cycles = re.findall(r"glwu\.grlc_2p5km\.t(\d{2})z\.grib2\.idx", s.get(f"{GLWU}/glwu.{day}/", timeout=60).text)
            if cycles:
                cycle = max(cycles)
                break
        else:
            return None
        url = f"{GLWU}/glwu.{day}/glwu.grlc_2p5km.t{cycle}z.grib2"
        idx = [line.split(":") for line in s.get(url + ".idx", timeout=60).text.splitlines()]
        step = lambda text: 0 if text == "anl" else int(text.split()[0]) if text.endswith("hour fcst") else None
        want = [k for k, e in enumerate(idx) if e[3] in GLWU_FIELDS and step(e[5]) is not None
                and (int(cycle) + step(e[5])) % CHART_STEP_H == 0]
        runs = []  # consecutive messages in one request
        for k in want:
            if runs and runs[-1][1] == k - 1:
                runs[-1][1] = k
            else:
                runs.append([k, k])
        fields, lat, lon = {}, None, None
        with tempfile.TemporaryFile() as f:
            for a, b in runs:
                end = int(idx[b + 1][1]) - 1 if b + 1 < len(idx) else ""
                r = s.get(url, headers={"Range": f"bytes={idx[a][1]}-{end}"}, timeout=120)
                r.raise_for_status()
                f.write(r.content)
            f.seek(0)
            for k in want:
                h = eccodes.codes_grib_new_from_file(f)
                try:
                    if lat is None:
                        lat, lon = eccodes.codes_get_array(h, "latitudes"), eccodes.codes_get_array(h, "longitudes")
                    v = eccodes.codes_get_values(h)
                    v = np.where(v == eccodes.codes_get(h, "missingValue"), np.nan, v)
                finally:
                    eccodes.codes_release(h)
                fields.setdefault(GLWU_FIELDS[idx[k][3]], {})[step(idx[k][5])] = v
    except Exception as e:
        print(f"  GLWU failed: {type(e).__name__}: {str(e)[:120]}", flush=True)
        return None
    hours = sorted(fields["hs"])
    water = np.isfinite(np.stack([fields["hs"][hr] for hr in hours])).any(axis=0)
    t0 = np.datetime64(f"{day[:4]}-{day[4:6]}-{day[6:]}T{cycle}:00")
    out = dict(times=np.array([t0 + np.timedelta64(hr, "h") for hr in hours], dtype="datetime64[m]"),
               lat=lat[water], lon=((lon[water] + 180) % 360) - 180)
    for name in GLWU_FIELDS.values():
        out[name] = np.stack([fields[name][hr][water] for hr in hours]) * (KT if name == "kt" else 1)
    print(f"  GLWU run {day} {cycle}z: {len(hours)} steps at {water.sum()} water points", flush=True)
    return out


def _glwu_point(g, lat, lon, reach=0.15):
    """Index of GLWU's nearest water point within reach degrees, or None."""
    d = np.hypot(g["lat"] - lat, (g["lon"] - lon) * np.cos(np.radians(lat)))
    k = int(np.argmin(d))
    return k if d[k] <= reach else None


def _glwu_rows(pts):
    """{key: rows} of GLWU's waves (one train: the whole sea, peak period) at the Great Lakes breaks."""
    from .points import lake_of

    g = _glwu()
    out = {}
    if g is None:
        return out
    hours = (g["times"][1] - g["times"][0]) / np.timedelta64(1, "h")
    for p in pts:
        if not lake_of(p["lat"], p["lon"]):
            continue
        k = _glwu_point(g, p["lat"], p["lon"])
        if k is None:
            continue
        hs, tp, d = g["hs"][:, k], g["tp"][:, k], g["dir"][:, k]
        nan = np.full(len(hs), np.nan)
        out[key(p)] = _rows(g["times"], p["lon"], [hs, tp, d, nan, nan, nan, hs, tp, d], hours)
    return out


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


def _from_chart(pts):
    """Waves and wind at each break from the swell chart (chart(), run first): ({key: wave rows},
    {key: wind rows}). All the waves stand in for the swell, from the nearest 1° cell with waves
    within 1.5°; the fallback when Open-Meteo is out of quota."""
    import base64
    import gzip

    m = json.loads(MAP_JSON.read_text())
    at, step = {n: k for k, (n, _) in enumerate(m["fields"])}, dict(m["fields"])
    cube = np.frombuffer(gzip.decompress(base64.b64decode(MAP_DATA.read_text())), np.uint8).reshape(len(m["times"]), len(at), m["nlat"], m["nlon"])
    times = np.array([np.datetime64(t.rstrip("Z")) for t in m["times"]])
    lat, lon = m["lat0"] + m["step"] * np.arange(m["nlat"]), m["lon0"] + m["step"] * np.arange(m["nlon"])
    hours = (times[1] - times[0]) / np.timedelta64(1, "h")
    wet = (cube[:, at["hs"]] != 255).any(axis=0)
    field = lambda n, a, b: np.where(cube[:, at[n], a, b] == 255, np.nan, cube[:, at[n], a, b] * step[n])
    waves, wind = {}, {}
    for p in pts:
        d = np.hypot(*np.meshgrid((lon - p["lon"] + 180) % 360 - 180, lat - p["lat"])[::-1] * np.array([1, np.cos(np.radians(p["lat"]))])[:, None, None])
        a, b = np.unravel_index(np.argmin(d), d.shape)
        nan = np.full(len(times), np.nan)
        air, rain = (_GFS[n][:, a, b] if n in _GFS else nan for n in ("air", "rain"))
        wind[key(p)] = _rows(times, p["lon"], [field("kt", a, b), field("wind_from", a, b), nan, air, rain, nan], hours)
        dw = np.where(wet, d, np.inf)
        a, b = np.unravel_index(np.argmin(dw), dw.shape)
        if dw[a, b] <= 1.5:
            hs, tp, dr = field("hs", a, b), field("tp", a, b), field("dir", a, b)
            waves[key(p)] = _rows(times, p["lon"], [hs, tp, dr, nan, nan, nan, hs, tp, dr], hours)
    return waves, wind
