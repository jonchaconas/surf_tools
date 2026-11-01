"""Regional pages: named coasts in Greece and on Lake Michigan.

Each coast is sampled from Open-Meteo at its own model cell (about 8 km), a few km
offshore. Enclosed seas rarely see groundswell, so besides the open-ocean good+ rule
these pages count "surfable windswell" mornings by a regional standard.
"""

import json
from functools import lru_cache

import numpy as np

from . import config, inland, rating

# A regional standard: total wave height min_hs-max_hs m with period >= min_period s, arriving
# within config.EXPOSURE degrees of the shore's facing. peak=True compares the dominant (peak)
# period, estimated from the mean wave period. Regions use WINDSWELL unless they set their own.
WINDSWELL = dict(noun="surfable", name="Surfable windswell", min_hs=0.9, max_hs=3.0, min_period=5.0, peak=False)
LAKE = dict(noun="good", name="Good lake surf", min_hs=1.22, max_hs=3.0, min_period=6.0, peak=True)  # 4 ft @ 6 s
# Open-ocean coasts use the good+ rule itself (rating.is_good on the swell, aimed at the shore).
OCEAN = dict(noun="good+", name="Good+ surf", ocean=True)

# extent (lon0, lon1, lat0, lat1); coasts are (name, area, lat, lon, facing), the point a few km offshore,
# optionally followed by a (lon, lat) for the map label where coasts crowd together.
REGIONS = {
    "Lake Michigan": dict(extent=(-88.4, -85.4, 41.45, 45.0),
                          standard=LAKE, coasts=[
        ("Whiting", "Indiana", 41.71, -87.48, 20, (-87.72, 41.58)),
        ("Gary", "Indiana", 41.68, -87.30, 0, (-87.38, 41.50)),
        ("Portage", "Indiana", 41.66, -87.17, 355, (-87.02, 41.50)),
        ("Chicago", "Illinois", 41.96, -87.60, 80, (-87.30, 42.12)),
        ("Racine", "Wisconsin", 42.74, -87.74, 90),
        ("Milwaukee", "Wisconsin", 43.06, -87.84, 90, (-88.05, 43.12)),
        ("Sheboygan", "Wisconsin", 43.75, -87.66, 90),
        ("New Buffalo", "Michigan", 41.83, -86.78, 320, (-86.50, 41.70)),
        ("St. Joseph", "Michigan", 42.11, -86.53, 280, (-86.30, 42.20)),
        ("Holland", "Michigan", 42.77, -86.25, 270),
        ("Grand Haven", "Michigan", 43.06, -86.29, 270),
        ("Ludington", "Michigan", 43.95, -86.50, 270),
        ("Frankfort", "Michigan", 44.63, -86.30, 270),
    ]),
    "Lake Superior": dict(extent=(-92.4, -84.3, 46.3, 49.1), standard=LAKE, coasts=[
        ("Park Point", "Duluth, Minnesota", 46.74, -92.04, 120),
        ("Stoney Point", "Minnesota", 46.91, -91.79, 135),
        ("Grand Marais MN", "Minnesota", 47.73, -90.33, 170),
        ("Little Girl's Point", "Michigan", 46.65, -90.32, 0),
        ("Ontonagon", "Michigan", 46.91, -89.33, 320),
        ("Marquette", "Michigan", 46.57, -87.36, 30),
        ("Au Train", "Michigan", 46.46, -86.84, 0),
        ("Grand Marais MI", "Michigan", 46.70, -85.98, 0, (-86.05, 46.42)),
        ("Whitefish Point", "Michigan", 46.79, -85.00, 330, (-85.25, 47.12)),
        ("Agawa Bay", "Ontario", 47.35, -84.72, 260),
        ("Terrace Bay", "Ontario", 48.75, -87.10, 190),
    ]),
    "Mediterranean": dict(extent=(-6.0, 36.5, 30.0, 46.0), coasts=[  # round the basin, Spain first
        ("Barcelona", "Spain", 41.37, 2.23, 130, (1.4, 40.75)),
        ("Marseille", "France", 43.25, 5.33, 230),
        ("Varazze", "Liguria, Italy", 44.33, 8.58, 180),
        ("Livorno", "Tuscany, Italy", 43.53, 10.26, 260),
        ("Santa Marinella", "Lazio, Italy", 42.02, 11.83, 230),
        ("Ajaccio", "Corsica, France", 41.90, 8.62, 250),
        ("Capo Mannu", "Sardinia, Italy", 40.03, 8.34, 280),
        ("Mondello", "Sicily, Italy", 38.22, 13.33, 20),
        ("Kastro", "Kyllini, Greece", 37.88, 21.03, 250),
        ("Kalo Nero", "Messinia, Greece", 37.30, 21.63, 260),
        ("Falasarna", "Crete, Greece", 35.50, 23.52, 270),
        ("Antalya", "Turkey", 36.84, 30.62, 190),
        ("Paphos", "Cyprus", 34.76, 32.38, 250),
        ("Beirut", "Lebanon", 33.89, 35.44, 270),
        ("Gaza", "Palestine", 31.52, 34.40, 295),
        ("Alexandria", "Egypt", 31.23, 29.90, 340),
        ("Susa", "Cyrenaica, Libya", 32.93, 21.97, 0),
        ("Benghazi", "Libya", 32.12, 20.00, 280),
        ("Tripoli", "Libya", 32.94, 13.18, 0),
        ("Hammamet", "Tunisia", 36.37, 10.66, 140),
        ("Kelibia", "Cap Bon, Tunisia", 36.85, 11.15, 80),
        ("Bizerte", "Tunisia", 37.30, 9.86, 15),
        ("Tabarka", "Tunisia", 36.99, 8.76, 0),
        ("Algiers", "Algeria", 36.80, 3.06, 0),
    ]),
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
    "Caspian Sea": dict(extent=(46.5, 55.0, 36.3, 47.2), coasts=[
        ("Aktau", "Kazakhstan", 43.63, 51.10, 240),
        ("Awaza", "Turkmenistan", 40.02, 52.90, 260),
        ("Babolsar", "Iran", 36.75, 52.66, 0),
        ("Ramsar", "Iran", 36.95, 50.65, 10),
        ("Lankaran", "Azerbaijan", 38.75, 48.93, 90),
        ("Shikhov", "Baku, Azerbaijan", 40.29, 49.83, 120),
        ("Novkhani", "Absheron, Azerbaijan", 40.60, 49.76, 340),
        ("Nabran", "Azerbaijan", 41.85, 48.72, 50),
        ("Derbent", "Dagestan, Russia", 42.06, 48.38, 80),
        ("Makhachkala", "Dagestan, Russia", 42.98, 47.58, 80),
    ]),
    "Black Sea": dict(extent=(27.0, 42.0, 40.6, 47.0), coasts=[
        ("Kilyos", "Istanbul, Turkey", 41.27, 29.03, 0),
        ("Şile", "Istanbul, Turkey", 41.19, 29.61, 15),
        ("Batumi", "Georgia", 41.65, 41.58, 290),
        ("Sochi", "Russia", 43.55, 39.68, 220),
        ("Anapa", "Russia", 44.89, 37.25, 270),
        ("Odesa", "Ukraine", 46.45, 30.82, 110),
        ("Vama Veche", "Romania", 43.75, 28.63, 95),
        ("Varna", "Bulgaria", 43.20, 28.00, 90),
        ("Sozopol", "Bulgaria", 42.43, 27.75, 70),
    ]),
    "Madagascar": dict(extent=(42.3, 51.6, -26.2, -11.4), standard=OCEAN, coasts=[
        ("Lavanono", "Deep south", -25.45, 44.95, 190),
        ("Fort Dauphin", "Tolagnaro, southeast", -25.05, 47.03, 120),
        ("Anakao", "Toliara, southwest", -23.66, 43.60, 250),
        ("Ifaty", "Toliara, southwest", -23.15, 43.55, 270),
        ("Andavadoaka", "Morombe, west", -22.07, 43.18, 270),
        ("Morondava", "West", -20.28, 44.20, 260),
        ("Maintirano", "West", -18.06, 43.96, 270),
        ("Mahajanga", "Northwest", -15.68, 46.26, 310),
        ("Nosy Be", "Northwest", -13.40, 48.18, 220),
        ("Ramena", "Antsiranana, north", -12.20, 49.42, 70),
        ("Sainte-Marie", "Nosy Boraha, east", -17.00, 49.92, 90),
        ("Mahambo", "East", -17.48, 49.52, 80),
        ("Manakara", "Southeast", -22.15, 48.07, 100),
    ]),
    "Peru": dict(extent=(-82.2, -70.0, -18.6, -3.3), standard=OCEAN, coasts=[  # north to south
        ("Máncora", "Piura", -4.08, -81.11, 290, (-80.45, -3.85)),
        ("Cabo Blanco", "Piura", -4.23, -81.28, 300, (-81.75, -4.0)),
        ("Lobitos", "Piura", -4.44, -81.33, 290, (-81.75, -4.75)),
        ("Pacasmayo", "La Libertad", -7.42, -79.62, 250, (-80.6, -7.15)),
        ("Chicama", "La Libertad", -7.72, -79.49, 250, (-80.6, -7.75)),
        ("Huanchaco", "Trujillo", -8.10, -79.17, 250, (-80.3, -8.35)),
        ("Huarmey", "Áncash", -10.10, -78.21, 243),
        ("Miraflores", "Lima", -12.17, -77.07, 220, (-78.4, -11.85)),
        ("Punta Hermosa", "Lima", -12.38, -76.86, 220, (-78.4, -12.45)),
        ("Cerro Azul", "Lima", -13.07, -76.51, 220, (-77.9, -13.25)),
        ("Camaná", "Arequipa", -16.66, -72.72, 193, (-74.2, -16.45)),
        ("Mollendo", "Arequipa", -17.06, -72.05, 220, (-73.6, -17.25)),
        ("Ilo", "Moquegua", -17.67, -71.39, 240, (-72.8, -17.95)),
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


@lru_cache(maxsize=None)
def _load(region, month):
    """Per coast, its morning arrays (points.mornings: MFWAM waves with both swell trains, ERA5
    wind) for one calendar month's seasons; None until every coast is downloaded."""
    from . import points

    out = []
    for _, _, lat, lon, *_ in REGIONS[region]["coasts"]:
        m = points.mornings(lat, lon)
        if m is None:
            return None
        out.append(points.season(m, (month,)))
    return out


def standard(region):
    return REGIONS[region].get("standard", WINDSWELL)


def _tally(m, facing, tp_ratio, std):
    """Mornings meeting the region's standard ("surf") and good+ at a break ("good": either swell
    train, light or offshore wind), and what the former looked like: total sea for the windswell
    and lake standards, which don't use wind, and the qualifying swell for the ocean one."""
    good = rating.break_good(m, facing)
    hs = m["hs"]
    if std.get("ocean"):
        surf, shown = good, rating.good_swell(m, facing)
    else:
        tm = m["tm"] * (tp_ratio if std["peak"] else 1.0)
        with np.errstate(invalid="ignore"):
            surf = (hs >= std["min_hs"]) & (hs <= std["max_hs"]) & (tm >= std["min_period"]) & \
                rating.exposed(m["dir"], facing)
        shown = (hs, tm, m["dir"])
    return dict(surf=int(surf.sum()), good=int(good.sum()), n=len(hs),
                hs=shown[0][surf], tm=shown[1][surf], dir=shown[2][surf])


def summary(region, tp_ratio):
    """Per coast: surfable % by month (Jan..Dec, NaN while downloading), year-round surfable and
    ocean good+ %, and the median surfable morning."""
    from . import points
    from .analysis import _circmean

    monthly = {m: _load(region, m) for m in range(1, 13)}
    rows = []
    for k, (name, area, lat, lon, facing, *label) in enumerate(REGIONS[region]["coasts"]):
        parts = {m: _tally(d[k], facing, tp_ratio, standard(region)) for m, d in monthly.items() if d is not None}
        surf, good, n = (sum(p[f] for p in parts.values()) for f in ("surf", "good", "n"))
        hs, tm, dr = (np.concatenate([p[f] for p in parts.values()]) if parts else np.array([])
                      for f in ("hs", "tm", "dir"))
        med = lambda a: float(np.median(a)) if len(a) else np.nan
        rows.append(dict(
            name=name, region=area, lat=lat, lon=lon, facing=facing, label=label[0] if label else None,
            months=[100 * parts[m]["surf"] / max(parts[m]["n"], 1) if m in parts else np.nan for m in range(1, 13)],
            pct=100 * surf / max(n, 1), ocean=100 * good / max(n, 1), complete=len(parts) == 12,
            hs=med(hs), tp=med(tm), dir=_circmean(dr),
            wind=points.typical_wind(points.season(m, range(1, 13))) if (m := points.mornings(lat, lon)) else None))
    return rows


def _km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, (a[0], a[1], b[0], b[1]))
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(h))


@lru_cache(maxsize=None)
def best(region, n, sep_km=160):
    """Names of the region's n best coasts year-round by its own standard, at least sep_km apart."""
    chosen = []
    for r in sorted(summary(region, inland.TP_RATIO), key=lambda r: -r["pct"]):
        if all(_km((r["lat"], r["lon"]), (c["lat"], c["lon"])) >= sep_km for c in chosen):
            chosen.append(r)
            if len(chosen) == n:
                break
    return tuple(c["name"] for c in chosen)


def coast_pct(region, months, std, tp_ratio):
    """{coast: (lat, lon, % of mornings meeting std)} over the given months; {} while any is downloading."""
    data = [_load(region, m) for m in months]
    if any(d is None for d in data):
        return {}
    out = {}
    for k, (name, _, lat, lon, facing, *_) in enumerate(REGIONS[region]["coasts"]):
        parts = [_tally(d[k], facing, tp_ratio, std) for d in data]
        out[name] = (lat, lon, 100 * sum(p["surf"] for p in parts) / max(sum(p["n"] for p in parts), 1))
    return out
