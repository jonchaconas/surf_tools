"""Turn morning samples into % good+ days: a 0.25 degree global field and per-spot figures."""

from functools import lru_cache

import cartopy.io.shapereader as shp
import numpy as np
from scipy.interpolate import RegularGridInterpolator, griddata
from scipy.ndimage import distance_transform_edt
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


def _pct(good, valid, axis):
    n = valid.sum(axis)
    pct = 100 * good.sum(axis) / np.maximum(n, 1)
    return np.where(n > 0.5 * n.max(), pct, np.nan)


def _upsample(pct, lat, lon):
    """Bilinear 1 degree -> 0.25 degree, plus where the coarse field had data."""
    has = np.isfinite(pct).astype(float)
    filled = _fill_nearest(pct)
    wrap = lambda a: np.concatenate([a[:, -1:], a, a[:, :1]], axis=1)
    lon_w = np.concatenate([lon[-1:] - 360, lon, lon[:1] + 360])
    la, lo = np.meshgrid(np.clip(LAT, lat.min(), lat.max()), LON, indexing="ij")
    f = RegularGridInterpolator((lat, lon_w), wrap(filled))((la, lo))
    h = RegularGridInterpolator((lat, lon_w), wrap(has))((la, lo))
    return f, h > 0.5


class Month:
    """Everything the PDF needs for one calendar month."""

    def __init__(self, month):
        self.month = month
        self.ocean = ww3.load(month)
        self.basins = inland.load(month)
        self.days = ww3.days(month)

    def field(self):
        """(% good+ on the 0.25 degree grid with land filled in, mask of water without data)."""
        ds = self.ocean
        good = rating.is_good(ds.shgt.values, ds.sper.values, ds.whgt.values)
        pct = _pct(good, np.isfinite(ds.shgt.values), axis=0)
        f, has = _upsample(pct, ds.lat.values, ds.lon.values)
        grid = np.where(has, f, np.nan)

        geoms = inland.geometries()
        for name, b in self.basins.items():
            p = self.basin_pct(name)
            ok = np.isfinite(p)
            la0, la1, lo0, lo1 = (lambda d, *bb: bb)(*inland.BASINS[name])
            rows = (LAT > la0) & (LAT < la1)
            cols = (LON > lo0) & (LON < lo1)
            la, lo = np.meshgrid(LAT[rows], LON[cols], indexing="ij")
            inside = contains_xy(geoms[name], lo, la)
            pts = (b["lat"][ok], b["lon"][ok])
            lin = griddata(pts, p[ok], (la[inside], lo[inside]), method="linear")
            near = griddata(pts, p[ok], (la[inside], lo[inside]), method="nearest")
            sub = grid[np.ix_(rows, cols)]
            sub[inside] = np.where(np.isnan(lin), near, lin)
            grid[np.ix_(rows, cols)] = sub

        nodata = water_mask() & np.isnan(grid)
        return _fill_nearest(grid), nodata

    def basin_pct(self, name):
        b = self.basins[name]
        return _pct(rating.is_good(b["hs"], b["tp"], b["ws"]), np.isfinite(b["hs"]), axis=1)

    def spots(self):
        """Per named break: % of mornings that are good+ AND arrive from the break's swell window."""
        geoms = inland.geometries()
        ds = self.ocean
        valid = (np.isfinite(ds.shgt.values).mean(0) > 0.5)
        la, lo = np.meshgrid(ds.lat.values, ds.lon.values, indexing="ij")
        ocean_pts = (la[valid], lo[valid], np.argwhere(valid))

        rows = []
        for name, region, lat, lon, facing in SPOTS:
            basin = next((k for k, g in geoms.items() if g.distance(Point(lon, lat)) < 1.0), None)
            if basin and basin in self.basins:
                b = self.basins[basin]
                k = _seaward(b["lat"], b["lon"], lat, lon, facing)
                if k is None:
                    continue
                hs, tp, dr, ws = b["hs"][k], b["tp"][k], b["dir"][k], b["ws"][k]
                source = "Open-Meteo"
            else:
                k = _seaward(ocean_pts[0], ocean_pts[1], lat, lon, facing)
                if k is None:
                    continue
                i, j = ocean_pts[2][k]
                pick = lambda v: ds[v].values[:, i, j]
                hs, tp, dr, ws = pick("shgt"), pick("sper"), pick("sdir"), pick("whgt")
                source = "WW3"
            good = rating.is_good(hs, tp, ws) & rating.exposed(dr, facing)
            n = np.isfinite(hs).sum()
            rows.append(dict(name=name, region=region, lat=lat, lon=lon, basin=basin,
                             pct=100 * good.sum() / max(n, 1), days=int(n), source=source))
        return sorted(rows, key=lambda r: -r["pct"])


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
