"""Export the binned morning data that the interactive atlas (web/) reads.

Each calendar month's mornings are counted, per grid cell and per named point, into bins of
height x period x cleanliness (and, for points, how far off the shore's facing the waves
arrive). The page turns any definition of "good" into sums over those bins, so people can
set their own criteria without the browser touching five years of daily model data.

grid.b64.txt and points.b64.txt hold gzip-compressed binary as base64 text (artifact hosting
serves text, not arbitrary binaries). Decompressed layout (little-endian):
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

from . import analysis, config, inland, rating, regions, tides, wind, ww3
from .spots import SPOTS
from .stretches import STRETCHES

OUT = config.ROOT / "web" / "data"
H_BINS, T_BINS = 20, 16               # 1 ft bins (0..19+), 1 s bins (<6 s, 6..20+ s)
CLEAN = [0.2, 0.35, 0.5, 0.75, 1.0]   # wind-sea / swell height edges -> 6 bins
ANGLE = [30, 45, 60, 75, 90]          # degrees off the shore's facing -> 6 bins
DEDUPE_KM = 40                        # a named break this close to a regional coast is drawn once


def _codes(hs, tp, ws=None, frm=None, facing=None):
    """Bin codes for morning samples. ws=None means total waves, where cleanliness doesn't apply."""
    h = np.clip(np.floor(np.nan_to_num(hs) * rating.FT), 0, H_BINS - 1).astype(np.int32)
    t = np.clip(np.floor(np.nan_to_num(tp)) - 5, 0, T_BINS - 1).astype(np.int32)
    if ws is None:
        c = np.zeros_like(h)
    else:
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.nan_to_num(ws, nan=0.0) / np.where(np.nan_to_num(hs) > 0, hs, np.nan)
        c = np.digitize(np.nan_to_num(ratio, nan=9.0), CLEAN)
    code = (h * T_BINS + t) * 6 + c
    if facing is not None:
        diff = np.abs((np.asarray(frm, dtype=float) - facing + 180) % 360 - 180)
        code = code * 6 + np.digitize(np.nan_to_num(diff, nan=180.0), ANGLE, right=True)
    return code


def _hist(codes, valid):
    u, n = np.unique(codes[valid], return_counts=True)
    return u.astype(np.uint16), n.astype(np.uint8), int(valid.sum())


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


def _grid(stride):
    """Ocean cells from WW3 (every stride-th 1 degree cell) plus the enclosed-sea points."""
    per_month, lat, lon = [], None, None
    for m in range(1, 13):
        ds = ww3.load(m).isel(lat=slice(None, None, stride), lon=slice(None, None, stride))
        hs, tp, ws = ds.shgt.values, ds.sper.values, ds.whgt.values
        lat, lon = ds.lat.values, ds.lon.values
        per_month.append((_codes(hs, tp, ws), rating.has_sea(hs, ws)))
        print(f"  grid month {m}", flush=True)
    days = sum(v.sum(0) for _, v in per_month)
    keep = np.argwhere(days > 0.5 * days.max())
    cells = [dict(lat=float(lat[i]), lon=float(lon[j]), size=float(stride), kind="ocean") for i, j in keep]
    hists = [[_hist(c[:, i, j], v[:, i, j]) for c, v in per_month] for i, j in keep]

    tp_ratio = inland.TP_RATIO
    basins = [inland.load(m) for m in range(1, 13)]
    for name, (spacing, *_) in inland.BASINS.items():
        if any(name not in b for b in basins):
            continue
        b0 = basins[0][name]
        for k in range(len(b0["lat"])):
            per = []
            for b in (bm[name] for bm in basins):
                hs = b["wave_hs"][k]
                per.append(_hist(_codes(hs, b["wave_tm"][k] * b["ratio"]), np.isfinite(hs)))
            if sum(p[2] for p in per):
                cells.append(dict(lat=float(b0["lat"][k]), lon=float(b0["lon"][k]), size=spacing, kind="enclosed"))
                hists.append(per)
    return cells, hists


def _points():
    """Named breaks, stretches and regional coasts, each with its own exposure-aware histograms."""
    coasts = [(region, c) for region, spec in regions.REGIONS.items() for c in spec["coasts"]]
    km = lambda a, b: regions._km(a, b)
    points, hists = [], []

    month_grids = []
    for m in range(1, 13):
        ds = ww3.load(m)
        details = {r: analysis._Grid(d) for r in ww3.DETAIL if (d := ww3.load(m, r)) is not None}
        month_grids.append((analysis._Grid(ds), details, inland.load(m)))

    for name, region, lat, lon, facing in SPOTS:
        if any(km((lat, lon), (c[2], c[3])) < DEDUPE_KM for _, c in coasts):
            continue  # the regional page's coast stands in for it
        basin = analysis._basin_of(lat, lon)
        per = []
        for globe, details, basins in month_grids:
            if basin:
                b = basins.get(basin)
                k = analysis._seaward(b["lat"], b["lon"], lat, lon, facing) if b else None
                if k is None:
                    per.append((np.zeros(0, np.uint16), np.zeros(0, np.uint8), 0))
                    continue
                hs, tp, frm = b["wave_hs"][k], b["wave_tm"][k] * b["ratio"], b["wave_dir"][k]
                per.append(_hist(_codes(hs, tp, None, frm, facing), np.isfinite(hs)))
            else:
                g = next((details[r] for r, (la0, la1, lo0, lo1) in ww3.DETAIL.items()
                          if r in details and la0 <= lat <= la1 and lo0 <= lon % 360 <= lo1), globe)
                k = analysis._seaward(g.lat, g.lon, lat, lon, facing)
                i, j = g.cells[k]
                hs, tp, frm, ws = g.hs[:, i, j], g.tp[:, i, j], g.frm[:, i, j], g.ws[:, i, j]
                per.append(_hist(_codes(hs, tp, ws, frm, facing), rating.has_sea(hs, ws)))
        points.append(dict(name=name, region=region, group="Named breaks", kind="break",
                           lat=lat, lon=lon, facing=facing, measure="waves" if basin else "swell"))
        hists.append(per)

    for num, (name, region, lat, lon, facing, (clat, clon)) in enumerate(STRETCHES, 1):
        per = []
        for globe, _, _ in month_grids:
            k = int(np.argmin(np.abs(globe.lat - clat) + np.abs(globe.lon - clon)))
            i, j = globe.cells[k]
            hs, tp, frm, ws = globe.hs[:, i, j], globe.tp[:, i, j], globe.frm[:, i, j], globe.ws[:, i, j]
            per.append(_hist(_codes(hs, tp, ws, frm, facing), rating.has_sea(hs, ws)))
        points.append(dict(name=name, region=region, group="Under-the-radar", kind="stretch", num=num,
                           lat=lat, lon=lon, facing=facing, measure="swell"))
        hists.append(per)

    for region, (name, area, lat, lon, facing, *_) in coasts:
        ocean = regions.standard(region).get("ocean")
        k = [c[0] for c in regions.REGIONS[region]["coasts"]].index(name)
        per = []
        for m in range(1, 13):
            data = regions._load(region, m)
            if data is None:
                per.append((np.zeros(0, np.uint16), np.zeros(0, np.uint8), 0))
                continue
            d = data[k]
            peak, mean = d["swell_wave_peak_period"], d["swell_wave_period"]
            if ocean:
                tp = np.where(np.isnan(peak) & np.isfinite(mean), mean * inland.TP_RATIO, peak)
                hs, frm, ws = d["swell_wave_height"], d["swell_wave_direction"], d["wind_wave_height"]
                per.append(_hist(_codes(hs, tp, ws, frm, facing), rating.has_sea(hs, ws)))
            else:
                hs = d["wave_height"]
                per.append(_hist(_codes(hs, d["wave_period"] * inland.TP_RATIO, None, d["wave_direction"], facing),
                                 np.isfinite(hs)))
        points.append(dict(name=name, region=area, group=region, kind="coast", lat=lat, lon=lon, facing=facing,
                           measure="swell" if ocean else "waves"))
        hists.append(per)

    for p in points:
        p["tide"] = tides.text(p["name"], rating.FT) or "none"
        p["wind"] = [list(np.round(w, 1)) if (w := wind.typical(p["name"], (m,))) else None for m in range(1, 13)]
    return points, hists


def _land():
    """Land outlines (Great Lakes cut out) as compact GeoJSON for drawing."""
    from .plot import _land as land_feature

    geom = next(iter(land_feature().geometries())).simplify(0.02, preserve_topology=True)
    gj = mapping(geom)

    def rnd(c):
        return [rnd(x) for x in c] if isinstance(c[0], (list, tuple)) else [round(c[0], 3), round(c[1], 3)]

    gj["coordinates"] = rnd(gj["coordinates"])
    return {"type": "Feature", "properties": {}, "geometry": gj}


def export(res=2):
    """Write web/data/{meta.json, grid.b64.txt, points.b64.txt, land.json}."""
    OUT.mkdir(parents=True, exist_ok=True)
    cells, grid = _grid(int(res))
    (OUT / "grid.b64.txt").write_text(_pack(grid, len(cells)))
    points, pts = _points()
    (OUT / "points.b64.txt").write_text(_pack(pts, len(points)))
    (OUT / "land.json").write_text(json.dumps(_land(), separators=(",", ":")))
    meta = dict(
        generated=date.today().isoformat(), res=res,
        seasons={m: list(config.seasons(m)) for m in range(1, 13)},
        bins=dict(height_ft=H_BINS, period_s=T_BINS, period_start=5, clean=CLEAN, angle=ANGLE),
        cells=dict(lat=[c["lat"] for c in cells], lon=[c["lon"] for c in cells],
                   size=[c["size"] for c in cells], enclosed=[c["kind"] == "enclosed" for c in cells]),
        points=points,
    )
    (OUT / "meta.json").write_text(json.dumps(meta, separators=(",", ":"), default=float))
    for f in sorted(OUT.iterdir()):
        print(f"  {f.name}: {f.stat().st_size / 1e6:.1f} MB")
