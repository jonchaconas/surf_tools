"""Export the data that the interactive atlas (web/) reads.

Grid cells (and the coastal band's source cells) have each calendar month's mornings counted into
bins of height x period x cleanliness, so the page turns any definition of "good" into sums over
bins without touching five years of daily model data. Named breaks, stretches and regional coasts
are few enough to ship their raw mornings instead (see _points), so every setting applies exactly.

The .b64.txt files hold gzip-compressed binary as base64 text (artifact hosting serves text, not
arbitrary binaries). Binned layout, decompressed (little-endian):
    Uint16 days[n * 12]     mornings with data per (entity, month)
    Uint16 entries[n * 12]  non-empty bins per (entity, month)
    Uint16 codes[total]     bin codes, entity-major then month
    Uint8  counts[total]    mornings in each bin
"""

import base64
import gzip
import json
from datetime import date

import numpy as np
from shapely.geometry import mapping

from . import analysis, config, inland, rating, regions, tides, ww3
from .spots import SPOTS
from .stretches import STRETCHES

OUT = config.ROOT / "web" / "data"
H_BINS, T_BINS = 20, 16               # 1 ft bins (0..19+), 1 s bins (<6 s, 6..20+ s)
CLEAN = [0.2, 0.35, 0.5, 0.75, 1.0]   # wind-sea / swell height edges -> 6 bins
DEDUPE_KM = 40                        # a named break this close to a regional coast is drawn once


def _codes(hs, tp, ws=None):
    """Bin codes for morning samples. ws=None means total waves, where cleanliness doesn't apply."""
    h = np.clip(np.floor(np.nan_to_num(hs) * rating.FT), 0, H_BINS - 1).astype(np.int32)
    t = np.clip(np.floor(np.nan_to_num(tp)) - 5, 0, T_BINS - 1).astype(np.int32)
    if ws is None:
        c = np.zeros_like(h)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.nan_to_num(ws, nan=0.0) / np.where(np.nan_to_num(hs) > 0, hs, np.nan)
        c = np.digitize(np.nan_to_num(ratio, nan=9.0), CLEAN)
    return (h * T_BINS + t) * 6 + c


def _hist(codes, valid):
    """(codes, counts, mornings). A count over 255 is split into repeated entries of one code."""
    u, n = np.unique(codes[valid], return_counts=True)
    us, ns = [u], [np.minimum(n, 255)]
    while (n > 255).any():
        u, n = u[n > 255], n[n > 255] - 255
        us.append(u)
        ns.append(np.minimum(n, 255))
    return np.concatenate(us).astype(np.uint16), np.concatenate(ns).astype(np.uint8), int(valid.sum())


def _pack(hists, n):
    """hists[entity][month] = (codes, counts, days) -> gzip bytes in the layout above."""
    days = np.zeros(n * 12, np.uint16)
    entries = np.zeros(n * 12, np.uint16)
    codes, counts = [], []
    for e in range(n):
        for m in range(12):
            c, k, d = hists[e][m]
            days[e * 12 + m], entries[e * 12 + m] = d, len(c)
            codes.append(c)
            counts.append(k)
    body = days.tobytes() + entries.tobytes() + np.concatenate(codes).tobytes() + np.concatenate(counts).tobytes()
    return base64.b64encode(gzip.compress(body, 6)).decode("ascii")


def _sector(frm):
    """Which of 16 compass sectors (0 = N) the waves come from; NaN -> 0 (those bins never pass)."""
    return (np.floor((np.nan_to_num(frm) % 360 + 11.25) / 22.5).astype(np.int32)) % 16


def _grid(stride):
    """Ocean cells from WW3 (every stride-th 1 degree cell) plus the enclosed-sea points.

    Returns (cells, samples, years): samples(e, m) gives cell e's mornings in calendar month m
    (0-11) as (codes, direction codes = code * 16 + sector, valid, year), from which _hists makes
    histograms over every year or a single one."""
    per_month, lat, lon = [], None, None
    for m in range(1, 13):
        ds = ww3.load(m).isel(lat=slice(None, None, stride), lon=slice(None, None, stride))
        hs, tp, ws = ds.shgt.values, ds.sper.values, ds.whgt.values
        lat, lon = ds.lat.values, ds.lon.values
        code = _codes(hs, tp, ws)
        per_month.append((code.astype(np.uint16), (code * 16 + _sector(ds.sdir.values)).astype(np.uint16),
                          rating.has_sea(hs, ws), ds.day.dt.year.values))
        print(f"  grid month {m}", flush=True)
    days = sum(v.sum(0) for _, _, v, _ in per_month)
    keep = np.argwhere(days > 0.1 * days.max())  # seasonal sea ice: rated over its open-water mornings
    cells = [dict(lat=float(lat[i]), lon=float(lon[j]), size=float(stride), kind="ocean", basin=None)
             for i, j in keep]
    where = [("ocean", int(i), int(j)) for i, j in keep]

    basins = [inland.load(m) for m in range(1, 13)]
    for name, (spacing, *_) in inland.BASINS.items():
        if any(name not in b for b in basins):
            continue
        b0 = basins[0][name]
        for k in range(len(b0["lat"])):
            per = []
            for b in (bm[name] for bm in basins):
                hs = b["wave_hs"][k]
                code = _codes(hs, b["wave_tm"][k] * b["ratio"])
                per.append((code, code * 16 + _sector(b["wave_dir"][k]), np.isfinite(hs),
                            b["day"].astype("datetime64[Y]").astype(int) + 1970))
            if sum(p[2].sum() for p in per):
                cells.append(dict(lat=float(b0["lat"][k]), lon=float(b0["lon"][k]), size=spacing, kind="enclosed",
                                  basin=name))
                where.append(("basin", per))

    def samples(e, m):
        w = where[e]
        if w[0] == "basin":
            return w[1][m]
        c, d, v, y = per_month[m]
        return c[:, w[1], w[2]], d[:, w[1], w[2]], v[:, w[1], w[2]], y

    years = sorted({int(y) for *_, ys in per_month for y in ys})
    return cells, samples, years


def _hists(samples, entities, year=None, direction=False):
    """Per entity, per calendar month: _hist of its codes (or direction codes), all years or one."""
    out = []
    for e in entities:
        per = []
        for m in range(12):
            c, d, v, y = samples(e, m)
            per.append(_hist(d if direction else c, v if year is None else v & (y == year)))
        out.append(per)
    return out


# One morning at a point, as bytes (255 = none): swell 1 and 2 as height (0.05 m), peak period
# (0.1 s) and direction (2 degrees); all waves as the same three; wind speed (0.5 kt) and direction;
# the year (since 2000); the tide, as % of the way from the break's low (5th percentile) to high (95th).
MORNING = [("hs1", 0.05), ("tp1", 0.1), ("dir1", 2), ("hs2", 0.05), ("tp2", 0.1), ("dir2", 2),
           ("hs", 0.05), ("tp", 0.1), ("dir", 2), ("kt", 0.5), ("wind_from", 2), ("year", 1), ("tide", 1)]


def _records(s, tide_range):
    """points.season() arrays -> uint8 records, one row per morning."""
    from .inland import TP_RATIO

    lo, hi = tide_range
    cols = []
    for key, step in MORNING:
        v = s["tm"] * TP_RATIO if key == "tp" else \
            s["day"].astype("datetime64[Y]").astype(int) + 1970 - 2000 if key == "year" else \
            np.clip(100 * (s["tide"] - lo) / max(hi - lo, 0.01), 0, 100) if key == "tide" else s[key]
        q = np.round(np.asarray(v, float) / step)
        cols.append(np.where(np.isfinite(q), np.clip(q, 0, 254), 255).astype(np.uint8))
    return np.stack(cols, axis=1)


def _forecast_extras(e):
    """What the page and forecast.py need for the coming week: the CDIP site (California) and the
    tide constituents {name: [amplitude m, phase deg]}, so the page can predict any morning's tide."""
    from . import cdip

    out = {}
    if (site := cdip.site_for(e["lat"], e["lon"])):
        out["cdip"] = site
    if (cons := tides._constituents(e["lat"], e["lon"])):
        out["tc"] = {c: [round(a, 4), round(ph, 2)] for c, (a, ph) in cons.items()}
    return out


def _point_records(points, hour):
    """The marked spots' records (as in _points) sampled at another local hour, same order and layout."""
    from . import points as P

    counts, recs = [], []
    for e in points:
        m = P.mornings(e["lat"], e["lon"], hour)
        for month in range(1, 13):
            r = _records(P.season(m, (month,)), e["tide_m"]) if m is not None else np.zeros((0, len(MORNING)), np.uint8)
            counts.append(len(r))
            recs.append(r)
    body = np.array(counts, np.uint16).tobytes() + np.concatenate(recs).tobytes()
    return base64.b64encode(gzip.compress(body, 6)).decode("ascii")


def _buoy_checks():
    """{(kind, name): its best buoy comparison (surfmap/buoys.py)}, for breaks and regional coasts:
    the pair with the most matched mornings among buoys within 40 km, else the nearest."""
    path = config.DATA / "buoys" / "validation.json"
    if not path.exists():
        return {}
    out, rows = {}, json.loads(path.read_text())
    for kind, name in {(v.get("kind", "break"), v["break"]) for v in rows}:
        pairs = [v for v in rows if v.get("kind", "break") == kind and v["break"] == name]
        close = [v for v in pairs if v["km"] <= 40]
        v = max(close, key=lambda v: v["n"]) if close else min(pairs, key=lambda v: v["km"])
        out[(kind, name)] = dict(buoy=v["buoy"], name=v["buoy_name"], km=round(v["km"]), n=v["n"], bias=round(v["bias_m"], 2),
                         r=round(v["r"], 2), trust=v["trust"], buoy4ft=round(v["pct4ft_buoy"]), model4ft=round(v["pct4ft_model"]),
                         nearshore=v["source"].startswith("CDIP"))
    return out


def _points():
    """Named breaks, stretches and regional coasts with their raw mornings (points.py), so the page
    applies any swell size, period, angle and wind cutoff exactly. Points still downloading are
    left out. Layout of points.b64.txt (gzip, base64): Uint16 mornings[n * 12] per point and
    month, then Uint8 records[total * len(MORNING)], point-major, month order."""
    from . import points as P

    coasts = [(region, c) for region, spec in regions.REGIONS.items() for c in spec["coasts"]]
    entries = []
    for name, region, lat, lon, facing in SPOTS:
        if any(regions._km((lat, lon), (c[2], c[3])) < DEDUPE_KM for _, c in coasts):
            continue  # the regional page's coast stands in for it
        basin = analysis._basin_of(lat, lon)
        entries.append(dict(name=name, region=region, group="Named breaks", kind="break", lat=lat, lon=lon,
                            facing=facing, measure="waves" if basin else "swell"))
    for num, (name, region, lat, lon, facing, _) in enumerate(STRETCHES, 1):
        entries.append(dict(name=name, region=region, group="Under-the-radar", kind="stretch", num=num,
                            lat=lat, lon=lon, facing=facing, measure="swell"))
    for region, (name, area, lat, lon, facing, *_) in coasts:
        entries.append(dict(name=name, region=area, group=region, kind="coast", lat=lat, lon=lon, facing=facing,
                            measure="swell" if regions.standard(region).get("ocean") else "waves"))

    from . import sst
    from .spots import WINDOWS

    checks = _buoy_checks()
    out, counts, recs = [], [], []
    for e in entries:
        m = P.mornings(e["lat"], e["lon"], config.LOCAL_HOUR)
        if m is None:
            print(f"  not downloaded yet, left out: {e['name']}", flush=True)
            continue
        tide = np.nanpercentile(m["tide"], [5, 95]) if np.isfinite(m["tide"]).any() else (0.0, 0.0)
        for month in range(1, 13):
            r = _records(P.season(m, (month,)), tide)
            counts.append(len(r))
            recs.append(r)
        e["tide_m"] = [round(float(x), 2) for x in tide]
        e.update(_forecast_extras(e))
        e["sst"] = sst.monthly(e["lat"], e["lon"])
        if e["kind"] == "break" and e["name"] in WINDOWS:
            e["window"] = list(WINDOWS[e["name"]])
        if (e["kind"], e["name"]) in checks:
            e["check"] = checks[(e["kind"], e["name"])]
        e["cell"] = [round(float(x), 3) for x in m["cell"]]
        e["source"] = m["source"]
        e["tide"] = tides.text(e["name"], rating.FT) or "none"
        out.append(e)
    body = np.array(counts, np.uint16).tobytes() + np.concatenate(recs).tobytes()
    return out, base64.b64encode(gzip.compress(body, 6)).decode("ascii")


GAP_KM = 100  # every stretch of shore gets a clickable spot within this distance: one per ~200 km of coast


def _coast_points(shore, existing):
    """[[name, region and country, lat, lon]]: points on the shoreline (the coastal band's pixels next
    to land) placed so no shore is more than GAP_KM from a spot, counting the named breaks,
    stretches, regional coasts and OpenStreetMap spots already there. Named after the nearest town."""
    from scipy.spatial import cKDTree

    from .discover import _xyz

    lat, lon = shore
    xyz = _xyz(lat, lon)
    have = cKDTree(_xyz(np.array([e[0] for e in existing]), np.array([e[1] for e in existing])))
    far = have.query(xyz)[0] * 6371
    chosen = np.empty((0, 3))
    picks = []
    for i in np.flatnonzero(far > GAP_KM)[np.argsort(-far[far > GAP_KM])]:  # biggest gaps first
        if len(chosen) and np.min(np.linalg.norm(chosen - xyz[i], axis=1)) * 6371 <= GAP_KM:
            continue
        chosen = np.vstack([chosen, xyz[i]])
        picks.append(i)
    places = _places()
    ptree = cKDTree(_xyz(np.array([p[3] for p in places]), np.array([p[4] for p in places])))
    out = []
    for i in picks:
        d, k = ptree.query(xyz[i])
        name, region, country = places[k][:3]
        out.append([name if d * 6371 <= 150 else f"{abs(lat[i]):.1f}°{'N' if lat[i] >= 0 else 'S'} "
                    f"{abs(lon[i]):.1f}°{'E' if lon[i] >= 0 else 'W'}",
                    ", ".join(x for x in (region, country) if x), round(float(lat[i]), 3), round(float(lon[i]), 3)])
    return out


TILE = 10  # degrees: detailed coastline tiles, fetched by the page only for the view when zoomed in


def _clockwise(geom):
    """Polygons (from any mix of geometries) with clockwise outer rings, as d3's spherical paths expect."""
    from shapely.geometry import MultiPolygon
    from shapely.geometry.polygon import orient

    parts = getattr(geom, "geoms", [geom])
    polys = [orient(g, sign=-1.0) for part in parts for g in getattr(part, "geoms", [part]) if g.geom_type == "Polygon"]
    return MultiPolygon(polys) if polys else MultiPolygon()


def _land_tiles():
    """Natural Earth 1:10m land (fill) and coastline (outline) cut into TILE-degree tiles, written to
    web/data/land10/{lat}_{lon}.json; returns the tile keys. Fills overlap their neighbours a little
    so no hairline shows between tiles; outlines come from the coastline set, so tile edges aren't
    drawn as coast."""
    import cartopy.io.shapereader as shp
    from shapely.geometry import box
    from shapely.ops import unary_union

    lakes = [r.geometry for r in shp.Reader(shp.natural_earth("10m", "physical", "lakes")).records()
             if r.attributes.get("name") in inland.GREAT_LAKES]
    land = unary_union(list(shp.Reader(shp.natural_earth("10m", "physical", "land")).geometries())).difference(unary_union(lakes))
    coast = unary_union(list(shp.Reader(shp.natural_earth("10m", "physical", "coastline")).geometries()))
    rnd = lambda c: [rnd(x) for x in c] if isinstance(c[0], (list, tuple)) else [round(c[0], 3), round(c[1], 3)]
    out = OUT / "land10"
    out.mkdir(exist_ok=True)
    for f in out.glob("*.json"):
        f.unlink()
    keys = []
    for la in range(-60, 80, TILE):
        for lo in range(-180, 180, TILE):
            fill = land.intersection(box(lo - 0.02, la - 0.02, lo + TILE + 0.02, la + TILE + 0.02))
            line = coast.intersection(box(lo, la, lo + TILE, la + TILE))
            if fill.is_empty and line.is_empty:
                continue
            fill = _clockwise(fill.simplify(0.004, preserve_topology=True))  # d3 reads anticlockwise rings as "all but this"
            gj = {k: mapping(g) for k, g in (("fill", fill), ("coast", line.simplify(0.004))) if not g.is_empty}
            for g in gj.values():
                g["coordinates"] = rnd(g["coordinates"]) if g.get("coordinates") else g.get("coordinates", [])
            (out / f"{la}_{lo}.json").write_text(json.dumps(gj, separators=(",", ":")))
            keys.append(f"{la}_{lo}")
    return keys


# Place names in Palestine (the atlas owner's choice): the country is Palestine, towns go by their
# Arabic names, and towns founded after 1948 by the Arab town they replaced or grew beside
# (Tel Aviv -> Yafa). The coast also gets the historic Arab towns along it, so coast points there
# are named after them.
PALESTINE = {"Haifa": ("Hayfa", 32.822, 34.978), "Nazareth": ("An-Nasira", 32.704, 35.296), "Tel Aviv": ("Yafa", 32.054, 34.752),
             "Ramla": ("Ar-Ramla", 31.927, 34.866), "Jerusalem": ("Al-Quds", 31.776, 35.234), "Beer Sheva": ("Bir as-Sab'", 31.244, 34.79)}
PALESTINE_COAST = [("'Akka", 32.927, 35.074), ("Al-Tantura", 32.607, 34.918), ("Qisarya", 32.497, 34.891), ("Umm Khalid", 32.330, 34.860),
                   ("Isdud", 31.754, 34.664), ("Al-Majdal 'Asqalan", 31.668, 34.574), ("Khan Yunis", 31.346, 34.306),
                   ("Rafah", 31.287, 34.245), ("Umm al-Rashrash", 29.557, 34.952)]


def _places():
    """[[name, region, country, lat, lon]] for every Natural Earth populated place (10m), so the page
    can name any stretch of coast after its nearest town (see PALESTINE)."""
    import cartopy.io.shapereader as shp

    out = []
    for r in shp.Reader(shp.natural_earth("10m", "cultural", "populated_places_simple")).records():
        a = r.attributes
        name, region, country, lat, lon = a["name"], a.get("adm1name") or "", a.get("adm0name") or "", a["latitude"], a["longitude"]
        if country in ("Israel", "Palestine", "Gaza Strip", "West Bank"):
            name, lat, lon = PALESTINE.get(name, (name, lat, lon))
            region, country = "", "Palestine"
        out.append([name, region, country, round(lat, 3), round(lon, 3)])
    return out + [[n, "", "Palestine", la, lo] for n, la, lo in PALESTINE_COAST]


def _land():
    """Land outlines (Great Lakes cut out) as compact GeoJSON for drawing."""
    from .plot import _land as land_feature

    geom = next(iter(land_feature().geometries())).simplify(0.02, preserve_topology=True)
    gj = mapping(geom)

    def rnd(c):
        return [rnd(x) for x in c] if isinstance(c[0], (list, tuple)) else [round(c[0], 3), round(c[1], 3)]

    gj["coordinates"] = rnd(gj["coordinates"])
    return {"type": "Feature", "properties": {}, "geometry": gj}


RASTER = analysis.RES   # 0.25 degrees: the page's smoothed field and coastal band live on this grid
COAST_DEG = 0.75        # water this close to land is "coast" and judged by its own exposure
RAY_KM, BASIN_RAY_KM, LAND_RES = 200, 100, 0.05   # open water a sector needs; enclosed seas build waves over shorter fetch
RAYS = 9               # sight lines across each 22.5 degree sector
ISLAND_KM2 = 1000      # land smaller than this only half-blocks a ray: swell wraps around it


def _land_raster():
    """Land at LAND_RES degrees for casting exposure rays (Great Lakes and Caspian as water):
    0 water, 1 an island under ISLAND_KM2, 2 larger land."""
    from PIL import Image, ImageDraw
    from scipy.ndimage import label

    from .plot import _land as land_feature

    w, h = int(360 / LAND_RES), int(180 / LAND_RES)
    img = Image.new("1", (w, h), 0)
    draw = ImageDraw.Draw(img)
    to_px = lambda ring: [((x + 180) / LAND_RES, (90 - y) / LAND_RES) for x, y in ring.coords]
    geom = next(iter(land_feature().geometries()))
    for poly in getattr(geom, "geoms", [geom]):
        draw.polygon(to_px(poly.exterior), fill=1)
        for hole in poly.interiors:
            draw.polygon(to_px(hole), fill=0)
    land = np.array(img, dtype=bool)
    parts, _ = label(land, structure=np.ones((3, 3)))
    px_km2 = (LAND_RES * 111.2) ** 2 * np.cos(np.radians(90 - (np.arange(h) + 0.5) * LAND_RES))
    area = np.bincount(parts.ravel(), weights=np.repeat(px_km2.astype(np.float32), w))
    area[0] = np.inf
    return np.where(land, np.where(area[parts] < ISLAND_KM2, 1, 2), 0).astype(np.uint8)


def _open_sectors(lat, lon, land, ray_km):
    """How open each point is to swell from each compass sector (k * 22.5 degrees, 0 = N), 0-255:
    the mean over RAYS sight lines across the sector, each 1 when it runs ray_km (per point) over
    water, 0 when larger land is in the way and 0.5 when only small islands are, since swell wraps
    around those. Bays and headlands close what they block; island chains leave windows."""
    steps = np.arange(5, RAY_KM + 1, 5)
    within = steps[None, :] <= np.asarray(ray_km)[:, None]
    h, w = land.shape
    out = np.zeros((len(lat), 16), np.float32)
    coslat = np.cos(np.radians(lat))[:, None]
    for k in range(16):
        for off in (np.arange(RAYS) + 0.5) / RAYS * 22.5 - 11.25:
            th = np.radians(k * 22.5 + off)
            la = lat[:, None] + steps[None, :] * np.cos(th) / 111.2
            lo = lon[:, None] + steps[None, :] * np.sin(th) / (111.2 * np.maximum(coslat, 0.05))
            r = np.clip(((90 - la) / LAND_RES).astype(int), 0, h - 1)
            c = ((lo + 180) / LAND_RES).astype(int) % w
            worst = (land[r, c] * within).max(axis=1)
            out[:, k] += np.where(worst == 2, 0.0, np.where(worst == 1, 0.5, 1.0))
    return np.round(out / RAYS * 255).astype(np.uint8)


def _raster(cells, res, land, charted):
    """Everything the page needs to draw a smooth field that follows the coast, on the 0.25 degree grid:
    cat: 0 ocean (interpolate the ocean grid), 1 enclosed sea (inverse-distance from its own points),
         2 water no source reaches (shown as no data), 3 land.
    The coastal band (water within COAST_DEG of land) overrides that with each pixel's own exposure:
    swell only counts from sectors open to it, read from its nearest seaward source cell. Its weight
    (alpha) fades from 1 at the shore toward the band's outer edge, so it meets the open-water field
    without a seam. A strip of land along the shore copies the nearest coastal value, so contours
    meet the coastline. Each coastal pixel also records how far it is from the nearest charted
    spot (named breaks, regional coasts, OpenStreetMap surf spots), for the uncharted highlight."""
    from scipy.ndimage import binary_dilation, distance_transform_edt
    from scipy.spatial import cKDTree
    from shapely import contains_xy

    from .discover import _bearing, _xyz

    LAT, LON = analysis.LAT, analysis.LON
    nlat, nlon = len(LAT), len(LON)
    la, lo = np.meshgrid(LAT, LON, indexing="ij")
    water = analysis.water_mask()

    ocean = [i for i, c in enumerate(cells) if c["kind"] == "ocean"]
    otree = cKDTree(_xyz(np.array([cells[i]["lat"] for i in ocean]), np.array([cells[i]["lon"] for i in ocean])))
    reach = 1.5 * res * 111.2
    near = otree.query(_xyz(la.ravel(), lo.ravel()))[0].reshape(nlat, nlon) * 6371 <= reach

    cat = np.where(water, 0, 3).astype(np.uint8)
    basin_of = np.full((nlat, nlon), -1, np.int16)
    names = [n for n in inland.BASINS if any(c["basin"] == n for c in cells)]
    bpix, bnb, bw = [], [], []
    for bi, name in enumerate(names):
        inside = binary_dilation(contains_xy(inland.geometries()[name], lo, la), iterations=4)
        basin_of[inside & (basin_of < 0)] = bi
        members = [i for i, c in enumerate(cells) if c["basin"] == name]
        ii = np.flatnonzero(inside.ravel())
        cl = np.cos(np.radians(np.mean([cells[i]["lat"] for i in members])))
        pts = np.array([[cells[i]["lat"], cells[i]["lon"] * cl] for i in members])
        k = min(4, len(members))
        d, j = cKDTree(pts).query(np.c_[la.ravel()[ii], lo.ravel()[ii] * cl], k=k)
        d, j = d.reshape(len(ii), k), j.reshape(len(ii), k)
        wgt = 1 / np.maximum(d, 0.05) ** 2
        wgt /= wgt.sum(1, keepdims=True)
        for p, jj, ww in zip(ii, j, wgt):
            bpix.append(int(p))
            bnb.extend([members[x] for x in jj] + [members[jj[0]]] * (4 - k))
            bw.extend([round(float(x), 4) for x in ww] + [0.0] * (4 - k))
        cat[inside & water] = 1
    cat[water & ~near & (basin_of < 0)] = 2

    # Coastal band: exposure per pixel, then its nearest seaward source.
    to_land = distance_transform_edt(water)  # in pixels
    coast = water & (to_land * RASTER <= COAST_DEG) & (cat != 2)
    ci = np.flatnonzero(coast.ravel())
    clat, clon = la.ravel()[ci], lo.ravel()[ci]
    in_basin = basin_of.ravel()[ci]
    print(f"  coastal pixels: {len(ci)}; casting exposure rays", flush=True)
    masks = _open_sectors(clat, clon, land, np.where(in_basin >= 0, BASIN_RAY_KM, RAY_KM))
    src = np.full(len(ci), -1, np.int32)
    dist, cand = otree.query(_xyz(clat, clon), k=8)
    for n in range(len(ci)):
        if in_basin[n] >= 0:
            members = [i for i, c in enumerate(cells) if c["basin"] == names[in_basin[n]]]
            ml = np.array([cells[i]["lat"] for i in members])
            mo = np.array([cells[i]["lon"] for i in members])
            order = np.argsort((ml - clat[n]) ** 2 + ((mo - clon[n]) * np.cos(np.radians(clat[n]))) ** 2)[:8]
            options = [(members[o], ml[o], mo[o]) for o in order]
        else:
            options = [(ocean[c], cells[ocean[c]]["lat"], cells[ocean[c]]["lon"])
                       for c, dd in zip(cand[n], dist[n]) if dd * 6371 <= reach]
        if not options:
            continue
        src[n] = options[0][0]
        for cell, olat, olon in options:  # prefer a source out through an open sector
            sector = int(((_bearing(clat[n], clon[n], olat, olon) + 11.25) % 360) // 22.5)
            if masks[n, sector] >= 128:
                src[n] = cell
                break
    band = COAST_DEG / RASTER
    alpha = np.clip((band + 0.5 - to_land.ravel()[ci]) / (band - 1), 0, 1)
    ok = src >= 0
    ci, masks, src, alpha = ci[ok], masks[ok], src[ok], alpha[ok]
    ctree = cKDTree(_xyz(np.array([c[0] for c in charted]), np.array([c[1] for c in charted])))
    charted_km = np.minimum(ctree.query(_xyz(clat[ok], clon[ok]))[0] * 6371, 65535).astype(np.uint16)
    sources = sorted(set(src.tolist()))
    sidx = {c: k for k, c in enumerate(sources)}

    # Land strip along the shore copies the nearest coastal pixel.
    is_coast = np.zeros(nlat * nlon, bool)
    is_coast[ci] = True
    d, (ri, rj) = distance_transform_edt(~is_coast.reshape(nlat, nlon), return_indices=True)
    fringe = np.flatnonzero(((cat == 3) & (d <= 2)).ravel())
    nearest = (ri * nlon + rj).ravel()[fringe]
    pos = np.full(nlat * nlon, -1, np.int64)
    pos[ci] = np.arange(len(ci))
    fringe_to = pos[nearest]
    keep = fringe_to >= 0
    fringe, fringe_to = fringe[keep], fringe_to[keep]

    coast_bin = (ci.astype(np.uint32).tobytes() + fringe.astype(np.uint32).tobytes()
                 + fringe_to.astype(np.uint32).tobytes()
                 + np.array([sidx[c] for c in src], np.uint16).tobytes() + charted_km.tobytes() + masks.tobytes()
                 + np.round(alpha * 255).astype(np.uint8).tobytes())
    print(f"  coastal band: {len(ci)} pixels, {len(fringe)} shore strip, {len(sources)} source cells", flush=True)
    return dict(
        res=RASTER, nlat=nlat, nlon=nlon,
        cat=base64.b64encode(gzip.compress(cat.tobytes(), 9)).decode("ascii"),
        basin=dict(pix=bpix, nb=bnb, w=bw),
        coast=dict(n=len(ci), fringe=len(fringe), sources=sources,
                   data=base64.b64encode(gzip.compress(coast_bin, 9)).decode("ascii")),
    ), sources, (clat[ok][alpha == 1], clon[ok][alpha == 1])


def export(res=2):
    """Write web/data/{meta.json, grid.b64.txt, points.b64.txt, coast.b64.txt, land.json, raster.json, osm.json,
    places.json}."""
    OUT.mkdir(parents=True, exist_ok=True)
    for old in [*OUT.glob("grid_*.b64.txt"), *OUT.glob("coast_*.b64.txt")]:
        old.unlink()
    cells, samples, years = _grid(int(res))
    (OUT / "grid.b64.txt").write_text(_pack(_hists(samples, range(len(cells))), len(cells)))
    points, pts = _points()
    (OUT / "points.b64.txt").write_text(pts)
    for old in OUT.glob("points_*.b64.txt"):
        old.unlink()
    for hour in config.HOURS:  # the other times of day, loaded by the page when chosen
        if hour != config.LOCAL_HOUR:
            (OUT / f"points_{hour:02d}.b64.txt").write_text(_point_records(points, hour))
            print(f"  spots at {hour:02d}:00", flush=True)
    (OUT / "land.json").write_text(json.dumps(_land(), separators=(",", ":")))
    (OUT / "places.json").write_text(json.dumps(_places(), separators=(",", ":"), ensure_ascii=False))
    tiles = _land_tiles()
    from . import osm

    land = _land_raster()
    from . import sst

    surf = osm.spots(land, LAND_RES)
    (OUT / "osm.json").write_text(json.dumps([[*o, sst.monthly(o[1], o[2])] for o in surf], separators=(",", ":"), ensure_ascii=False))
    charted = [(la, lo) for _, _, la, lo, _ in SPOTS] + [(la, lo) for _, la, lo in surf] + \
        [(c[2], c[3]) for spec in regions.REGIONS.values() for c in spec["coasts"]]
    raster, sources, shore = _raster(cells, res, land, charted)
    gaps = [[*g, sst.monthly(g[2], g[3])] for g in _coast_points(shore, charted + [(e[2], e[3]) for e in STRETCHES])]
    (OUT / "coastpts.json").write_text(json.dumps(gaps, separators=(",", ":"), ensure_ascii=False))
    print(f"  coast points filling the gaps: {len(gaps)}", flush=True)
    (OUT / "raster.json").write_text(json.dumps(raster, separators=(",", ":")))
    (OUT / "coast.b64.txt").write_text(_pack(_hists(samples, sources, direction=True), len(sources)))
    # The same split by year, which the page loads only when the years slider narrows the record.
    for y in years:
        (OUT / f"grid_{y}.b64.txt").write_text(_pack(_hists(samples, range(len(cells)), y), len(cells)))
        (OUT / f"coast_{y}.b64.txt").write_text(_pack(_hists(samples, sources, y, True), len(sources)))
        print(f"  year {y}", flush=True)
    meta = dict(
        generated=date.today().isoformat(), res=res,
        seasons={m: list(config.seasons(m)) for m in range(1, 13)}, years=years, tiles=tiles, tile=TILE,
        hour=config.LOCAL_HOUR, hours=list(config.HOURS),
        bins=dict(height_ft=H_BINS, period_s=T_BINS, period_start=5, clean=CLEAN),
        morning=[k for k, _ in MORNING], scale=[step for _, step in MORNING],
        offshore=config.OFFSHORE,
        cells=dict(lat=[c["lat"] for c in cells], lon=[c["lon"] for c in cells],
                   size=[c["size"] for c in cells], enclosed=[c["kind"] == "enclosed" for c in cells]),
        ocean_res=res,
        points=points,
    )
    (OUT / "meta.json").write_text(json.dumps(meta, separators=(",", ":"), default=float))
    for f in sorted(OUT.iterdir()):
        print(f"  {f.name}: {f.stat().st_size / 1e6:.1f} MB")
