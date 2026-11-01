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

from . import config, inland, rating, ww3
from .spots import SPOTS

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
        basins={k: dict(lat=b["lat"], lon=b["lon"], n=rating.has_sea(b["hs"], b["ws"]).sum(1),
                        good=rating.is_good(b["hs"], b["tp"], b["ws"]).sum(1))
                for k, b in basins.items()},
        detail={r: g.counts(dds) for r, (dds, g) in details.items()},
        tp_ratio=next((b["tp_ratio"] for b in basins.values() if b["tp_ratio"]), None),
        days=ww3.days(month),
    )


def _spot_counts(globe, details, basins):
    """Good+ mornings with swell inside each break's exposure window, from its seaward cell.

    Breaks in an enclosed sea use Open-Meteo, breaks inside a DETAIL box its 0.5 degree
    grid, and the rest the global 1 degree grid.
    """
    out = {}
    for name, _, slat, slon, facing in SPOTS:
        basin = _basin_of(slat, slon)
        if basin:
            if basin not in basins:  # this month's Open-Meteo data not downloaded yet
                continue
            b = basins[basin]
            k = _seaward(b["lat"], b["lon"], slat, slon, facing)
            if k is None:
                continue
            h, t, d, w = b["hs"][k], b["tp"][k], b["dir"][k], b["ws"][k]
        else:
            g = next((details[r] for r, (la0, la1, lo0, lo1) in ww3.DETAIL.items()
                      if r in details and la0 <= slat <= la1 and lo0 <= slon % 360 <= lo1), globe)
            k = _seaward(g.lat, g.lon, slat, slon, facing)
            if k is None:
                continue
            i, j = g.cells[k]
            h, t, d, w = g.hs[:, i, j], g.tp[:, i, j], g.frm[:, i, j], g.ws[:, i, j]
        good = rating.is_good(h, t, w) & rating.exposed(d, facing)
        out[name] = (int(good.sum()), int(rating.has_sea(h, w).sum()))
    return out


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
        return _fill_nearest(grid), nodata

    def spots(self):
        """Per named break: % of mornings that are good+ AND arrive from the break's swell window."""
        rows = []
        for name, region, lat, lon, _ in SPOTS:
            basin = _basin_of(lat, lon)
            per = [c.spots.get(name) for c in self.counts]
            if (basin and basin not in self.ready) or None in per:
                continue
            good, n = sum(p[0] for p in per), sum(p[1] for p in per)
            rows.append(dict(name=name, region=region, lat=lat, lon=lon, basin=basin,
                             pct=100 * good / max(n, 1), days=n))
        return sorted(rows, key=lambda r: -r["pct"])

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
