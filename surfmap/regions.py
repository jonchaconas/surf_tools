"""Regional pages: named coasts in Greece and on Lake Michigan.

Each coast is sampled from Open-Meteo at its own model cell (about 8 km), a few km
offshore. Enclosed seas rarely see groundswell, so besides the open-ocean good+ rule
these pages count "surfable windswell" mornings by a regional standard.
"""

import json

import numpy as np

from . import config, inland, rating, ww3

# Surfable windswell: total wave height MIN_HS-MAX_HS m with mean period >= MIN_PERIOD s,
# arriving within config.EXPOSURE degrees of the shore's facing.
MIN_HS, MAX_HS, MIN_PERIOD = 0.9, 3.0, 5.0

# extent (lon0, lon1, lat0, lat1); coasts are (name, area, lat, lon, facing), the point a few km offshore.
REGIONS = {
    "Greece": dict(extent=(19.0, 29.6, 34.6, 41.0), coasts=[
        ("Agios Gordios", "Corfu, west coast", 39.53, 19.76, 260),
        ("Zakynthos west", "Ionian Islands", 37.82, 20.56, 250),
        ("Kastro", "Kyllini, NW Peloponnese", 37.88, 21.03, 250),
        ("Kalo Nero", "Messinia, W Peloponnese", 37.30, 21.63, 260),
        ("Falasarna", "West Crete", 35.50, 23.52, 270),
        ("Chania", "North Crete", 35.56, 23.97, 0),
        ("Kokkini Hani", "North Crete", 35.36, 25.25, 10),
        ("Komos", "South Crete", 34.99, 24.72, 225),
        ("Tinos north", "Cyclades", 37.66, 25.12, 340),
        ("Kymi", "East Evia", 38.62, 24.18, 70),
        ("Ialysos", "Rhodes, northwest", 36.43, 28.11, 300),
    ]),
    "Lake Michigan": dict(extent=(-88.4, -85.4, 41.45, 45.0), coasts=[
        ("Whiting", "Indiana", 41.71, -87.48, 20),
        ("Gary", "Indiana", 41.68, -87.30, 0),
        ("Portage", "Indiana", 41.66, -87.17, 355),
        ("Chicago", "Illinois", 41.96, -87.60, 80),
        ("Racine", "Wisconsin", 42.74, -87.74, 90),
        ("Milwaukee", "Wisconsin", 43.06, -87.84, 90),
        ("Sheboygan", "Wisconsin", 43.75, -87.66, 90),
        ("New Buffalo", "Michigan", 41.83, -86.78, 320),
        ("St. Joseph", "Michigan", 42.11, -86.53, 280),
        ("Holland", "Michigan", 42.77, -86.25, 270),
        ("Grand Haven", "Michigan", 43.06, -86.29, 270),
        ("Ludington", "Michigan", 43.95, -86.50, 270),
        ("Frankfort", "Michigan", 44.63, -86.30, 270),
    ]),
}
OUT = config.DATA / "regions"


def _path(region, year, month):
    return OUT / f"{region.lower().replace(' ', '_')}_{year}-{month:02d}.json"


def fetch():
    """One Open-Meteo request per region and month-season, under the same quota throttle as inland.py."""
    OUT.mkdir(parents=True, exist_ok=True)
    throttle, total = inland._Throttle(inland.PER_MINUTE, inland.PER_HOUR), 0.0
    for region, spec in REGIONS.items():
        lats, lons = [c[2] for c in spec["coasts"]], [c[3] for c in spec["coasts"]]
        for month in config.MONTHS:
            for year in config.seasons(month):
                path = _path(region, year, month)
                if path.exists():
                    continue
                d0, d1 = inland._window(year, month)
                weight = len(lats) * max(1.0, ((d1 - d0).days + 1) / 14)
                throttle.wait(weight)
                rows = inland._request(lats, lons, d0, d1)
                for row, a, b in zip(rows, lats, lons):
                    row["requested"] = [a, b]
                path.write_text(json.dumps(rows))
                total += weight
                print(f"  {path.name} ({total:.0f} weighted calls so far)", flush=True)


def _load(region, month):
    """Per coast, {variable: morning samples} for one calendar month; None until every season is cached."""
    paths = [_path(region, y, month) for y in config.seasons(month)]
    if not all(p.exists() for p in paths):
        return None
    keep = np.array(ww3.days(month), dtype="datetime64[D]")
    per = [{} for _ in REGIONS[region]["coasts"]]
    for p in paths:
        for k, row in enumerate(json.loads(p.read_text())):
            for v, a in inland.morning(row, keep).items():
                per[k].setdefault(v, []).append(a)
    return [{v: np.concatenate(a) for v, a in d.items()} for d in per]


def _tally(m, facing, tp_ratio):
    hs, tm, frm = m["wave_height"], m["wave_period"], m["wave_direction"]
    with np.errstate(invalid="ignore"):
        surf = (hs >= MIN_HS) & (hs <= MAX_HS) & (tm >= MIN_PERIOD) & rating.exposed(frm, facing)
    tp = m["swell_wave_peak_period"]
    if np.isnan(tp).all():  # no peak period on the Great Lakes, as in inland.load
        tp = m["swell_wave_period"] * tp_ratio
    good = rating.is_good(m["swell_wave_height"], tp, m["wind_wave_height"]) & \
        rating.exposed(m["swell_wave_direction"], facing)
    return dict(surf=int(surf.sum()), good=int(good.sum()), n=int(np.isfinite(hs).sum()),
                hs=hs[surf], tm=tm[surf], dir=frm[surf])


def summary(region, tp_ratio):
    """Per coast: surfable % by month (Jan..Dec, NaN while downloading), year-round surfable and
    ocean good+ %, and the median surfable morning."""
    from . import wind
    from .analysis import _circmean

    monthly = {m: _load(region, m) for m in range(1, 13)}
    rows = []
    for k, (name, area, lat, lon, facing) in enumerate(REGIONS[region]["coasts"]):
        parts = {m: _tally(d[k], facing, tp_ratio) for m, d in monthly.items() if d is not None}
        surf, good, n = (sum(p[f] for p in parts.values()) for f in ("surf", "good", "n"))
        hs, tm, dr = (np.concatenate([p[f] for p in parts.values()]) if parts else np.array([])
                      for f in ("hs", "tm", "dir"))
        med = lambda a: float(np.median(a)) if len(a) else np.nan
        rows.append(dict(
            name=name, region=area, lat=lat, lon=lon, facing=facing,
            months=[100 * parts[m]["surf"] / max(parts[m]["n"], 1) if m in parts else np.nan for m in range(1, 13)],
            pct=100 * surf / max(n, 1), ocean=100 * good / max(n, 1), complete=len(parts) == 12,
            hs=med(hs), tp=med(tm), dir=_circmean(dr), wind=wind.typical(name, tuple(range(1, 13)))))
    return rows
