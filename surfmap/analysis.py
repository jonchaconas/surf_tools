"""Turn morning samples into % good+ days: a 0.25 degree global field and per-spot figures.

Each calendar month is reduced once to counts (good mornings, valid mornings, summed
swell vectors). A Period adds up the counts of its months, so a season or full-year page
costs no extra loading and gives every morning equal weight.
"""

from dataclasses import dataclass
from functools import lru_cache
import calendar

import cartopy.io.shapereader as shp
import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.ndimage import binary_dilation, distance_transform_edt
from shapely import contains_xy
from shapely.geometry import Point
from shapely.ops import unary_union

from . import config, inland, points, rating, regions, ww3
from .spots import SPOTS, WINDOWS
from .stretches import STRETCHES

RES = 0.25
LAT = np.arange(-90 + RES / 2, 90, RES)
LON = np.arange(-180 + RES / 2, 180, RES)
MAX_SPOT_KM = 250


def _fill_nearest(a):
    """Replace NaNs with the nearest finite value, so contours run cleanly under the land mask."""
    idx = distance_transform_edt(np.isnan(a), return_distances=False, return_indices=True)
    return a[tuple(idx)]


@lru_cache
def water_mask():
    """True for 0.25 degree cells in the sea or the Great Lakes (Natural Earth 1:50m)."""
    path = config.DATA / "water_mask_025.npy"
    if path.exists():
        return np.load(path)
    ocean = unary_union(list(shp.Reader(shp.natural_earth("50m", "physical", "ocean")).geometries()))
    water = ocean.union(inland.geometries()["great_lakes"])
    la, lo = np.meshgrid(LAT, LON, indexing="ij")
    mask = contains_xy(water, lo, la)
    np.save(path, mask)
    return mask


@lru_cache
def _basin_of(lat, lon):
    """Which Open-Meteo basin a coastal point belongs to, if any."""
    return next((k for k, g in inland.geometries().items() if g.distance(Point(lon, lat)) < 1.0), None)


def _pct(good, n):
    """% good+, blank where a cell has under half the usual number of mornings (sea ice)."""
    return np.where(n > 0.5 * n.max(), 100 * good / np.maximum(n, 1), np.nan)


def _upsample(pct, lat, lon, rows=slice(None), cols=slice(None), wrap=True):
    """Bilinear coarse -> 0.25 degree over LAT[rows] x LON[cols], plus where the coarse field had data."""
    has = np.isfinite(pct).astype(float)
    filled = _fill_nearest(pct)
    if wrap:
        lon = np.concatenate([lon[-1:] - 360, lon, lon[:1] + 360])
        filled, has = (np.concatenate([a[:, -1:], a, a[:, :1]], axis=1) for a in (filled, has))
    la, lo = np.meshgrid(np.clip(LAT[rows], lat.min(), lat.max()),
                         np.clip(LON[cols], lon.min(), lon.max()), indexing="ij")
    f = RegularGridInterpolator((lat, lon), filled)((la, lo))
    h = RegularGridInterpolator((lat, lon), has)((la, lo))
    return f, h > 0.5


def _box_index(box):
    """Boolean row and column masks of the 0.25 degree grid inside a ww3.DETAIL box."""
    la0, la1, lo0, lo1 = box
    lo0, lo1 = (lo0 + 180) % 360 - 180, (lo1 + 180) % 360 - 180  # boxes don't cross the dateline
    return (LAT >= la0) & (LAT <= la1), (LON >= lo0) & (LON <= lo1)


def _seaward(plat, plon, lat, lon, facing):
    """Index of the nearest point within MAX_SPOT_KM, preferring points out to sea from the break."""
    dlon = (plon - lon + 180) % 360 - 180
    x, y = dlon * np.cos(np.radians(lat)), plat - lat
    km = 111.2 * np.hypot(x, y)
    bearing = np.degrees(np.arctan2(x, y)) % 360
    off = np.abs((bearing - facing + 180) % 360 - 180)
    for allowed in (np.where(off <= 90, km, np.inf), km):
        k = int(np.argmin(allowed))
        if allowed[k] <= MAX_SPOT_KM:
            return k
    return None


@dataclass
class Counts:
    """One calendar month across its seasons, reduced to sums."""
    lat: np.ndarray       # WW3 1 degree axes
    lon: np.ndarray
    good: np.ndarray      # good+ mornings per ocean cell
    n: np.ndarray         # mornings with data per ocean cell
    u: np.ndarray         # wave-power-weighted swell heading components
    v: np.ndarray
    w: np.ndarray         # total weight
    spots: dict           # spot name -> (good, n)
    stretches: dict       # stretch name -> dict(good, n) plus breaker height, period, direction of good mornings
    basins: dict          # complete Open-Meteo basins -> dict(lat, lon, good, n) per point
    detail: dict          # downloaded ww3.DETAIL boxes -> dict(lat, lon, good, n) at 0.5 degrees
    tp_ratio: float | None
    days: list


class _Grid:
    """A WW3 grid's morning samples, with its ocean cells listed for nearest-cell lookups."""

    def __init__(self, ds):
        self.hs, self.tp, self.frm, self.ws = (ds[v].values for v in ("shgt", "sper", "sdir", "whgt"))
        ocean = rating.has_sea(self.hs, self.ws).mean(0) > 0.5
        la, lo = np.meshgrid(ds.lat.values, ds.lon.values, indexing="ij")
        self.lat, self.lon, self.cells = la[ocean], lo[ocean], np.argwhere(ocean)

    def counts(self, ds):
        return dict(lat=ds.lat.values, lon=ds.lon.values, n=rating.has_sea(self.hs, self.ws).sum(0),
                    good=rating.is_good(self.hs, self.tp, self.ws).sum(0))


@lru_cache(maxsize=None)
def month_counts(month):
    ds = ww3.load(month)
    globe = _Grid(ds)
    w = np.nan_to_num(globe.hs**2 * globe.tp)
    heading = np.radians(globe.frm + 180)
    details = {}
    for region in ww3.DETAIL:
        if (dds := ww3.load(month, region)) is not None:
            details[region] = (dds, _Grid(dds))
    basins = inland.load(month)
    return Counts(
        lat=ds.lat.values, lon=ds.lon.values,
        good=globe.counts(ds)["good"], n=globe.counts(ds)["n"],
        u=np.nansum(w * np.sin(heading), 0), v=np.nansum(w * np.cos(heading), 0), w=w.sum(0),
        spots=_spot_counts(globe, {r: g for r, (_, g) in details.items()}, basins),
        stretches=_stretch_counts(globe),
        basins={k: dict(lat=b["lat"], lon=b["lon"], n=np.isfinite(b["wave_hs"]).sum(1),
                        good=_relaxed(b["wave_hs"], b["wave_tm"], b["ratio"]).sum(1))
                for k, b in basins.items()},
        detail={r: g.counts(dds) for r, (dds, g) in details.items()},
        tp_ratio=next((b["tp_ratio"] for b in basins.values() if b["tp_ratio"]), None),
        days=ww3.days(month),
    )


def _spot_counts(globe, details, basins):
    """Good mornings with waves inside each break's exposure window, from its seaward cell.

    Breaks in an enclosed sea use that basin's points and the relaxed lake standard; breaks
    inside a DETAIL box use its 0.5 degree grid; the rest use the global 1 degree grid and good+.
    """
    out = {}
    for name, _, slat, slon, facing in SPOTS:
        basin = _basin_of(slat, slon)
        if basin:
            if basin in basins:  # absent while this month's enclosed-sea data is still downloading
                b = basins[basin]
                k = _seaward(b["lat"], b["lon"], slat, slon, facing)
                if k is not None:
                    out[name] = _relaxed_tally(b, k, facing)
            continue
        g = next((details[r] for r, (la0, la1, lo0, lo1) in ww3.DETAIL.items()
                  if r in details and la0 <= slat <= la1 and lo0 <= slon % 360 <= lo1), globe)
        k = _seaward(g.lat, g.lon, slat, slon, facing)
        if k is None:
            continue
        i, j = g.cells[k]
        out[name] = _tally(g.hs[:, i, j], g.tp[:, i, j], g.frm[:, i, j], g.ws[:, i, j], facing)
    return out


def _stretch_counts(globe):
    """Each discovered stretch, read from the exact 1 degree cell discover.py chose for it."""
    out = {}
    for name, _, _, _, facing, (clat, clon) in STRETCHES:
        k = int(np.argmin(np.abs(globe.lat - clat) + np.abs(globe.lon - clon)))
        if abs(globe.lat[k] - clat) + abs(globe.lon[k] - clon) > 0.01:  # cell blank this month
            continue
        i, j = globe.cells[k]
        out[name] = _tally(globe.hs[:, i, j], globe.tp[:, i, j], globe.frm[:, i, j], globe.ws[:, i, j], facing)
    return out


def _relaxed(hs, tm, ratio):
    """Enclosed seas rarely see groundswell, so their contours use the regional lake standard
    (regions.LAKE: total waves 4-10 ft at a dominant period of at least 6 s) instead of good+."""
    std = regions.LAKE
    with np.errstate(invalid="ignore"):
        return (hs >= std["min_hs"]) & (hs <= std["max_hs"]) & (tm * (ratio if std["peak"] else 1) >= std["min_period"])


def _relaxed_tally(b, k, facing):
    """_tally's counterpart for a break in an enclosed sea, by the relaxed standard."""
    h, t, d = b["wave_hs"][k], b["wave_tm"][k] * b["ratio"], b["wave_dir"][k]
    good = _relaxed(b["wave_hs"][k], b["wave_tm"][k], b["ratio"]) & rating.exposed(d, facing)
    return dict(good=int(good.sum()), n=int(np.isfinite(h).sum()), hs=h[good],
                hb=rating.breaker_height(h[good], t[good]), tp=t[good], dir=d[good])


def _tally(h, t, d, w, facing):
    """Good+ mornings with swell inside the facing's exposure window, and what those mornings looked like."""
    good = rating.is_good(h, t, w) & rating.exposed(d, facing)
    return dict(good=int(good.sum()), n=int(rating.has_sea(h, w).sum()), hs=h[good],
                hb=rating.breaker_height(h[good], t[good]), tp=t[good], dir=d[good])


def _merge(per):
    """Add up _tally results over months: % good+ and the median good morning."""
    good, n = sum(p["good"] for p in per), sum(p["n"] for p in per)
    hs, hb, tp, dr = (np.concatenate([p[k] for p in per]) for k in ("hs", "hb", "tp", "dir"))
    med = lambda a: float(np.median(a)) if len(a) else np.nan
    return dict(pct=100 * good / max(n, 1), days=n, hs=med(hs), hb=med(hb), tp=med(tp), dir=_circmean(dr))


def _circmean(deg):
    if not len(deg):
        return np.nan
    r = np.radians(deg)
    return np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean())) % 360


class Period:
    """Everything one PDF page needs, for a single month or a run of months."""

    def __init__(self, months, title, table_title):
        self.months, self.title, self.table_title = tuple(months), title, table_title
        self.counts = [month_counts(m) for m in self.months]
        self.lat, self.lon = self.counts[0].lat, self.counts[0].lon
        self.ready = set.intersection(*(set(c.basins) for c in self.counts))
        self.detail = set.intersection(*(set(c.detail) for c in self.counts))
        self.tp_ratio = next((c.tp_ratio for c in self.counts if c.tp_ratio), None)
        self.days = sorted(d for c in self.counts for d in c.days)

    def _sum(self, attr):
        return sum(getattr(c, attr) for c in self.counts)

    def describe(self):
        if len(self.months) == 1:
            m = self.months[0]
            years = list(config.seasons(m))
            full = sum(calendar.monthrange(y, m)[1] for y in years)
            gap = "" if len(self.days) == full else f"; {full - len(self.days)} days missing from the archive"
            return f"{len(years)} seasons, {years[0]}–{years[-1]} ({len(self.days)} mornings{gap})"
        span = f"{self.days[0]:%b %Y} to {self.days[-1]:%b %Y}"
        return f"{len(self.days):,} mornings, {span}" + ("" if len(self.months) == 12 else ", these months only")

    def field(self):
        """(% good+ on the 0.25 degree grid with land filled in, mask of water without data)."""
        f, has = _upsample(_pct(self._sum("good"), self._sum("n")), self.lat, self.lon)
        for region in self.detail:
            parts = [c.detail[region] for c in self.counts]
            pct = _pct(sum(p["good"] for p in parts), sum(p["n"] for p in parts))
            rows, cols = _box_index(ww3.DETAIL[region])
            fr, hr = _upsample(pct, parts[0]["lat"], parts[0]["lon"], rows, cols, wrap=False)
            box = np.ix_(rows, cols)
            f[box] = np.where(hr, fr, f[box])
            has[box] |= hr
        grid = np.where(has, f, np.nan)

        geoms = inland.geometries()
        for name in self.ready:
            parts = [c.basins[name] for c in self.counts]
            p = _pct(sum(b["good"] for b in parts), sum(b["n"] for b in parts))
            ok = np.isfinite(p)
            _, la0, la1, lo0, lo1 = inland.BASINS[name]
            rows = (LAT > la0) & (LAT < la1)
            cols = (LON > lo0) & (LON < lo1)
            la, lo = np.meshgrid(LAT[rows], LON[cols], indexing="ij")
            inside = contains_xy(geoms[name], lo, la)
            pts = (parts[0]["lat"][ok], parts[0]["lon"][ok])
            lin = griddata(pts, p[ok], (la[inside], lo[inside]), method="linear")
            near = griddata(pts, p[ok], (la[inside], lo[inside]), method="nearest")
            sub = grid[np.ix_(rows, cols)]
            sub[inside] = np.where(np.isnan(lin), near, lin)
            grid[np.ix_(rows, cols)] = sub

        water = water_mask()
        nodata = water & np.isnan(grid)
        # Grow the mask onto land too, so the land layer's finer coastline leaves no slivers of
        # extrapolated colour along shores without data. Real water data is never covered.
        nodata |= binary_dilation(nodata, iterations=2) & ~water
        # Clamp: interpolation can leave -1e-19 and the like, which falls outside the 0% contour band.
        return np.clip(_fill_nearest(grid), 0, 100), nodata

    def spots(self):
        """Per named break: % of mornings good+ at the break itself (points.py: MFWAM at its ~8 km
        cell, either swell train inside its window, ERA5 morning wind); enclosed-sea breaks by
        the relaxed standard. Breaks still downloading are left out."""
        rows = []
        for name, region, lat, lon, facing in SPOTS:
            basin = _basin_of(lat, lon)
            s = points.summary(lat, lon, facing, self.months, relaxed=bool(basin), window=WINDOWS.get(name))
            if s is not None:
                rows.append(dict(name=name, region=region, lat=lat, lon=lon, basin=basin, **s))
        return sorted(rows, key=lambda r: -r["pct"])

    def stretches(self):
        """Per discovered stretch, numbered as in stretches.py: % good+ at its coast and its median good morning."""
        rows = []
        for num, (name, region, lat, lon, facing, _) in enumerate(STRETCHES, 1):
            s = points.summary(lat, lon, facing, self.months)
            if s is not None:
                rows.append(dict(num=num, name=name, region=region, lat=lat, lon=lon, **s))
        return rows

    def swell_arrows(self, step=10, min_consistency=0.3):
        """Predominant swell heading on a coarse open-ocean grid: (lon, lat, u, v) unit vectors.

        Mornings are summed as vectors pointing where the swell travels, weighted by wave
        power (Hs^2 * Tp). Cells whose mean vector is under min_consistency of the total
        weight have no predominant direction and are left out.
        """
        u, v, w, n = self._sum("u"), self._sum("v"), self._sum("w"), self._sum("n")
        with np.errstate(invalid="ignore", divide="ignore"):
            consistency = np.hypot(u, v) / w
        rows = np.isclose((self.lat + step / 4) % step, 0)  # -72.5, -62.5, ... 77.5
        cols = np.isclose(self.lon % step, step / 2)         # -175, -165, ... 175
        keep = np.zeros(n.shape, bool)
        keep[np.ix_(rows, cols)] = True
        keep &= (n > 0.5 * n.max()) & (consistency >= min_consistency)
        la, lo = np.meshgrid(self.lat, self.lon, indexing="ij")
        norm = np.hypot(u, v)[keep]
        return lo[keep], la[keep], u[keep] / norm, v[keep] / norm
