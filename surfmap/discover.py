"""Find under-the-radar surf coast in the WW3 record.

A stretch is a coastal model cell (within COAST_KM of land) whose mornings are often good+
with swell travelling toward that coast: MIN_YEAR % year-round, or MIN_PEAK % over its
best three consecutive months. It must be at least FAR_KM from every break in spots.py
(the well-known ones) and within TOWN_KM of a town of TOWN_POP people or more. Of stretches
closer than SEP_KM, only the best is kept.

`python -m surfmap discover` rewrites stretches.py from the current data.
"""

import cartopy.io.shapereader as shp
import numpy as np
from scipy.spatial import cKDTree

from . import config, rating, ww3
from .spots import SPOTS

COAST_KM, FAR_KM, SEP_KM, TOWN_KM, TOWN_POP = 120, 300, 200, 150, 5000
NAME_POP, NAME_COAST_KM = 1000, 30  # stretches are named after the nearest town on the coast
# Where the nearest seaside town in Natural Earth is across water or far off: generated -> (name, region).
RENAME = {
    "Tumby Bay coast": ("Kangaroo Island coast", "Australia"),
    "Pisco coast": ("Marcona coast", "Peru"),
    "Puerto Villamil coast": ("Isabela coast", "Galápagos, Ecuador"),
    "Puerto Baquerizo Moreno coast": ("San Cristóbal coast", "Galápagos, Ecuador"),
}
MIN_YEAR, MIN_PEAK = 40, 65
R = 6371.0


def _xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]


def _bearing(lat1, lon1, lat2, lon2):
    """Initial compass bearing from point 1 to point 2, degrees."""
    p1, p2, dl = np.radians(lat1), np.radians(lat2), np.radians(lon2 - lon1)
    y = np.sin(dl) * np.cos(p2)
    x = np.cos(p1) * np.sin(p2) - np.sin(p1) * np.cos(p2) * np.cos(dl)
    return np.degrees(np.arctan2(y, x)) % 360


def _coastline():
    pts = []
    for g in shp.Reader(shp.natural_earth("50m", "physical", "land")).geometries():
        for poly in getattr(g, "geoms", [g]):
            c = np.asarray(poly.exterior.coords)
            pts.append(c[:: max(1, len(c) // 4000)])
    return np.vstack(pts)


def _town_tree(towns):
    return cKDTree(_xyz(np.array([t["latitude"] for t in towns]), np.array([t["longitude"] for t in towns])))


def _coastal_cells(coast):
    """WW3 ocean cells near land, with the nearest coast point and the bearing out to sea."""
    ds = ww3.load(1)
    la, lo = np.meshgrid(ds.lat.values, ds.lon.values, indexing="ij")
    ocean = rating.has_sea(ds.shgt.values, ds.whgt.values).mean(0) > 0.5
    d, idx = cKDTree(_xyz(coast[:, 1], coast[:, 0])).query(_xyz(la[ocean], lo[ocean]))
    near = d * R <= COAST_KM
    rows, cols = (a[near] for a in np.nonzero(ocean))
    lat, lon, cp = la[ocean][near], lo[ocean][near], coast[idx[near]]
    return dict(lat=lat, lon=lon, rows=rows, cols=cols, coast_lat=cp[:, 1], coast_lon=cp[:, 0],
                seaward=_bearing(cp[:, 1], cp[:, 0], lat, lon))


def _monthly(cells):
    """(12, cell) % of mornings good+ with swell from the seaward side, plus year-round %."""
    pct, good, n = np.zeros((12, len(cells["lat"]))), 0, 0
    for m in range(1, 13):
        o = ww3.load(m)
        hs, tp, frm, ws = (o[v].values[:, cells["rows"], cells["cols"]] for v in ("shgt", "sper", "sdir", "whgt"))
        g = rating.is_good(hs, tp, ws) & rating.exposed(frm, cells["seaward"])
        k = rating.has_sea(hs, ws).sum(0)
        pct[m - 1] = 100 * g.sum(0) / np.maximum(k, 1)
        good, n = good + g.sum(0), n + k
        print(f"  month {m}", flush=True)
    return pct, 100 * good / np.maximum(n, 1)


def find():
    coast = _coastline()
    cells = _coastal_cells(coast)
    monthly, year = _monthly(cells)
    peak = np.max([np.roll(monthly, -s, 0)[:3].mean(0) for s in range(12)], axis=0)

    known = cKDTree(_xyz(np.array([s[2] for s in SPOTS]), np.array([s[3] for s in SPOTS])))
    far = known.query(_xyz(cells["lat"], cells["lon"]))[0] * R >= FAR_KM

    places = [r.attributes for r in shp.Reader(shp.natural_earth("10m", "cultural", "populated_places_simple")).records()]
    big = [t for t in places if (t.get("pop_max") or 0) >= TOWN_POP]
    tkm = _town_tree(big).query(_xyz(cells["coast_lat"], cells["coast_lon"]))[0] * R
    small = [t for t in places if (t.get("pop_max") or 0) >= NAME_POP]
    shore = cKDTree(_xyz(coast[:, 1], coast[:, 0])).query(
        _xyz(np.array([t["latitude"] for t in small]), np.array([t["longitude"] for t in small])))[0] * R
    seaside = [t for t, d in zip(small, shore) if d <= NAME_COAST_KM]
    seaside_tree = _town_tree(seaside)

    ok = ((year >= MIN_YEAR) | (peak >= MIN_PEAK)) & far & (tkm <= TOWN_KM)
    chosen = []
    for i in sorted(np.flatnonzero(ok), key=lambda i: -year[i]):
        p = _xyz(cells["lat"][i:i + 1], cells["lon"][i:i + 1])
        if all(np.linalg.norm(p - _xyz(cells["lat"][j:j + 1], cells["lon"][j:j + 1])) * R >= SEP_KM for j in chosen):
            chosen.append(i)

    town, used = {}, set()  # nearest seaside town not already naming a better stretch
    for i in chosen:
        near = seaside_tree.query(_xyz(cells["coast_lat"][i:i + 1], cells["coast_lon"][i:i + 1]), k=10)[1][0]
        town[i] = next((seaside[j] for j in near if seaside[j]["name"] not in used), seaside[near[0]])
        used.add(town[i]["name"])
    return [dict(name=f"{town[i]['name']} coast", region=town[i]["adm0name"],
                 lat=round(float(cells["coast_lat"][i]), 2), lon=round(float(cells["coast_lon"][i]), 2),
                 facing=int(round(cells["seaward"][i])), cell=(float(cells["lat"][i]), float(cells["lon"][i])),
                 year=float(year[i]), peak=float(peak[i]), town_km=float(tkm[i]),
                 monthly=[float(v) for v in monthly[:, i]])
            for i in chosen]


def write(stretches, path=config.ROOT / "surfmap" / "stretches.py"):
    for s in stretches:
        s["name"], s["region"] = RENAME.get(s["name"], (s["name"], s["region"]))
    lines = [
        '"""Under-the-radar surf coast found by discover.py (regenerate with `python -m surfmap discover`).',
        "",
        "Each entry is (name, region, lat, lon, facing, cell): a coast point, the compass direction",
        "it faces out to sea, and the WW3 1 degree cell its figures come from. Numbering on the",
        "maps follows this order (best year-round first).",
        '"""',
        "",
        "STRETCHES = [",
    ]
    for s in stretches:
        lines.append(f'    ({s["name"]!r}, {s["region"]!r}, {s["lat"]}, {s["lon"]}, {s["facing"]}, {s["cell"]}),'
                     f'  # {s["year"]:.0f}% year-round, {s["peak"]:.0f}% peak season')
    lines.append("]")
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {len(stretches)} stretches to {path}")
