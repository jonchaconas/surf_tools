"""Refraction and shoaling from each point's offshore wave cell in to the break, by backward ray tracing.

Breaks are rated from Copernicus MFWAM's ~8 km cell some way out to sea, so a break tucked into a bay
or round a point is rated on swell that never reaches it as such. This module estimates the change
between that cell and TARGET_DEPTH (10 m) off the break as a transfer table per point: for each
offshore swell direction (every 5 degrees) and peak period (6-20 s), the height coefficient
K = Kr * Ks and the direction the swell then comes from. In California CDIP's MOP model (cdip.py)
does this properly, and validate() compares the two there. Bathymetry is NOAA NCEI's ETOPO 2022
(15 arc-second surface elevation, over THREDDS OPeNDAP, no login), cached per point in data/bathy/;
tables go to data/refract/. `build()` makes both; `transfer(lat, lon)` applies a table.

Method: backward ray tracing (Dorrestein 1960; O'Reilly & Guza 1993, as CDIP's MOP does). From the
point nearest the break that is TARGET_DEPTH deep, a fan of rays (one every FAN degrees, all round)
is traced out to sea through the smoothed depth field, with linear dispersion for the phase speed
and Snell's law (d heading / ds from the speed gradient across the ray). A ray ends once it is as
far from its start as the model cell and as deep as it (or deep for its period), where the cell's
waves apply; it then carries the offshore direction that feeds it. A ray that runs aground or loops
is blocked. Along a ray E(f, theta) * c * cg is conserved, so the nearshore spectrum is the offshore
one read along the rays. For an offshore swell from D with a Gaussian spread(T),
    K^2 = sum over rays of  spread(theta_offshore - D) * FAN * (c cg)_offshore / (c cg)_TARGET_DEPTH,
which is Kr^2 Ks^2 for straight parallel contours (checked to ~1%), stays finite at caustics, and
goes to 0 where every ray from that direction is blocked by land (a headland's shadow). The
nearshore direction is the weighted mean of the rays' launch directions.

Against CDIP at the nine California breaks, the refracted direction is far closer than the
offshore one at eight (not Ocean Beach). The height fixes the sheltered breaks (Steamer Lane, Ventura, Malibu),
but the gains at open beaches (shoaling, focusing on long-period swell) overshoot. So use it with
max_k=1: the losses, not the gains.

Assumptions and limitations: linear waves at the peak period only, no currents, diffraction, bottom
friction or breaking. The offshore swell is taken as uniform where the rays end, so islands beyond
that are left to the model. ETOPO's 15 arc-second (~460 m) grid, smoothed by another cell, blurs
reefs, canyons and the 10 m contour; it is best on surveyed coasts (as in California) and often
interpolated elsewhere, occasionally walling off a bay (transfer() then returns None). The
directional spread is assumed (narrower for longer periods), not read from the model.
"""

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

import numpy as np
import requests
import xarray as xr
from scipy.ndimage import gaussian_filter

from . import config, points

ETOPO = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/15s/15s_surface_elev_netcdf/"
         "ETOPO_2022_v1_15s_{}_surface.nc")
PER_DEG = 240     # ETOPO 2022 15 arc-second cells per degree
TILE = 15         # degrees per ETOPO tile, named by its north-west corner
PAD = 0.4         # degrees of bathymetry each side of a point (longitude widened by 1/cos lat)
BATHY = config.DATA / "bathy"
OUT = config.DATA / "refract"
M_PER_DEG = 111_195.0
G = 9.81

TARGET_DEPTH = 10.0              # m, where the nearshore waves are given
DIRS = np.arange(0, 360, 5.0)    # offshore swell direction (from), degrees
PERIODS = np.arange(6, 21, 2.0)  # peak period, s
FAN = 0.25        # degrees between rays leaving the break
STEP = 100.0      # m along each ray (Runge-Kutta midpoint)
SMOOTH = 1.0      # Gaussian smoothing of the depth grid, in cells: rays need a smooth speed gradient
MIN_REACH = 3000  # m: rays run at least this far even when the model cell is nearer
MAX_REACH = 35000  # m: and stop by this far (inside the bathymetry box) however shallow it still is
DRY = 0.5         # m: a ray in water shallower than this has run aground
# A named break that rays say gets almost nothing from its whole swell window (12 s swell within
# config.EXPOSURE of its facing) means a bathymetry artefact, such as a bay ETOPO has walled off.
MIN_WINDOW_K = 0.3

# Lakes' surface heights above sea level, since ETOPO gives their beds: (lat0, lat1, lon0, lon1, m).
LEVELS = [(41.3, 42.95, -83.6, -78.8, 174.3),  # Lake Erie
          (43.15, 44.3, -79.9, -76.0, 74.8),   # Lake Ontario
          (46.3, 49.1, -92.4, -84.3, 183.4),   # Lake Superior
          (41.5, 46.2, -88.1, -79.7, 176.4),   # Lakes Michigan and Huron
          (36.5, 47.5, 46.5, 55.5, -27.5)]     # Caspian Sea


def spread(tp):
    """Assumed directional spread of a swell train (Gaussian s.d., degrees): narrower for longer periods."""
    return np.clip(30.0 - np.asarray(tp, float), 12.0, 24.0)


# --- Bathymetry -----------------------------------------------------------------------------

def _tile_name(row, col):
    """ETOPO tile holding global tile row/column (counted from 90S, 180W), e.g. N45W135 for 30-45N, 135-120W."""
    top, left = -90 + TILE * (row + 1), -180 + TILE * col
    return f"{'N' if top >= 0 else 'S'}{abs(top):02d}{'E' if left >= 0 else 'W'}{abs(left):03d}"


def _block(name, r0, r1, c0, c1):
    """Elevations [r0..r1, c0..c1] (rows south to north) of one tile, from OPeNDAP's text form."""
    url = ETOPO.format(name) + f".ascii?z[{r0}:1:{r1}][{c0}:1:{c1}]"
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=180)
            r.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    rows = re.findall(r"^\[\d+\], (.*)$", r.text.split("\n---", 1)[1], re.M)
    z = np.array([np.fromstring(row, sep=",") for row in rows], np.float32)
    if z.shape != (r1 - r0 + 1, c1 - c0 + 1):
        raise ValueError(f"ETOPO {name}: got {z.shape}")
    return z


def bathymetry(lat, lon):
    """(lats, lons, elevation in m) round a point: ETOPO 2022 15 arc-second surface elevation,
    cached in data/bathy/."""
    path = BATHY / f"{points._stem(lat, lon)}.npz"
    if not path.exists():
        n, last = TILE * PER_DEG, 360 * PER_DEG - 1
        dlon = PAD / np.cos(np.radians(lat))
        r0, r1 = (int((v + 90) * PER_DEG) for v in (lat - PAD, lat + PAD))
        c0, c1 = (min(max(int((v + 180) * PER_DEG), 0), last) for v in (lon - dlon, lon + dlon))
        z = np.empty((r1 - r0 + 1, c1 - c0 + 1), np.float32)
        for tr in range(r0 // n, r1 // n + 1):
            for tc in range(c0 // n, c1 // n + 1):
                a, b = max(r0, tr * n), min(r1, tr * n + n - 1)
                c, d = max(c0, tc * n), min(c1, tc * n + n - 1)
                z[a - r0:b - r0 + 1, c - c0:d - c0 + 1] = _block(_tile_name(tr, tc), a - tr * n, b - tr * n,
                                                                 c - tc * n, d - tc * n)
        BATHY.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".part.npz")
        np.savez_compressed(tmp, lat=-90 + (np.arange(r0, r1 + 1) + 0.5) / PER_DEG,
                            lon=-180 + (np.arange(c0, c1 + 1) + 0.5) / PER_DEG, z=z)
        tmp.rename(path)
    with np.load(path) as f:
        return f["lat"], f["lon"], f["z"]


def _level(lat, lon):
    return next((m for la0, la1, lo0, lo1, m in LEVELS if la0 <= lat <= la1 and lo0 <= lon <= lo1), 0.0)


def _grid(lat, lon):
    """Smoothed water depth h[y, x] (m, negative on land) on an east/north metre grid centred on
    the point: (x, y, h)."""
    lats, lons, z = bathymetry(lat, lon)
    h = np.maximum(_level(lat, lon) - z.astype(float), -5.0)  # else smoothing lets a cliff shoal the sea below it
    return (lons - lon) * M_PER_DEG * np.cos(np.radians(lat)), (lats - lat) * M_PER_DEG, gaussian_filter(h, SMOOTH)


def _xy(lat, lon, lat0, lon0):
    return (lon - lon0) * M_PER_DEG * np.cos(np.radians(lat0)), (lat - lat0) * M_PER_DEG


# --- Linear waves ---------------------------------------------------------------------------

def _k(h, T):
    """Wavenumber from omega^2 = g k tanh(k h), by Newton's method from k0 / sqrt(tanh(k0 h))."""
    w2 = (2 * np.pi / T) ** 2
    h = np.maximum(h, DRY)
    k = w2 / G / np.sqrt(np.tanh(w2 / G * h))
    for _ in range(4):
        t = np.tanh(k * h)
        k = k - (G * k * t - w2) / (G * t + G * k * h * (1 - t * t))
    return k


def _speeds(h, T):
    """Phase and group speed (m/s) at depth h for period T."""
    k = _k(h, T)
    kh = k * np.maximum(h, DRY)
    c = 2 * np.pi / T / k
    return c, c * (0.5 + kh / np.sinh(np.minimum(2 * kh, 600.0)))


# --- Rays -----------------------------------------------------------------------------------

def _start(x, y, h, px, py, facing, depth):
    """The point nearest (px, py) that is `depth` deep, looking first on the seaward side (within
    90 degrees of facing): the nearest deep-enough cell, then bisection back towards (px, py)."""
    X, Y = np.meshgrid(x, y)
    if _depth_at(x, y, h, px, py) >= depth:
        return px, py
    deep = h >= depth
    km = np.hypot(X - px, Y - py)
    off = np.abs((np.degrees(np.arctan2(X - px, Y - py)) - facing + 180) % 360 - 180)
    pick = deep & (off <= 90) if (deep & (off <= 90) & (km < 10_000)).any() else deep
    i = np.argmin(np.where(pick, km, np.inf))
    ex, ey = X.flat[i], Y.flat[i]
    lo, hi = 0.0, 1.0
    for _ in range(25):
        mid = (lo + hi) / 2
        if _depth_at(x, y, h, px + mid * (ex - px), py + mid * (ey - py)) >= depth:
            hi = mid
        else:
            lo = mid
    return px + hi * (ex - px), py + hi * (ey - py)


def _depth_at(x, y, h, px, py):
    fx, fy = np.interp(px, x, np.arange(len(x))), np.interp(py, y, np.arange(len(y)))
    i, j = min(int(fy), len(y) - 2), min(int(fx), len(x) - 2)
    wy, wx = fy - i, fx - j
    return float(h[i, j] * (1 - wy) * (1 - wx) + h[i, j + 1] * (1 - wy) * wx
                 + h[i + 1, j] * wy * (1 - wx) + h[i + 1, j + 1] * wy * wx)


def _fields(h, dx, dy):
    """Per period [c, dc/dx, dc/dy, h] on the grid, shape (periods, ny, nx, 4)."""
    out = np.empty((len(PERIODS), *h.shape, 4), np.float32)
    for p, T in enumerate(PERIODS):
        c = _speeds(h, T)[0]
        cy, cx = np.gradient(c, dy, dx)
        out[p] = np.stack([c, cx, cy, h], -1)
    return out


def _sample(F, p, px, py, x0, y0, dx, dy):
    """Bilinear [c, c_x, c_y, h] at each ray's position (p is its period index); inside=False off the grid."""
    fx, fy = (px - x0) / dx, (py - y0) / dy
    inside = (fx >= 0) & (fy >= 0) & (fx <= F.shape[2] - 1) & (fy <= F.shape[1] - 1)
    j = np.clip(fx.astype(int), 0, F.shape[2] - 2)
    i = np.clip(fy.astype(int), 0, F.shape[1] - 2)
    wx, wy = np.clip(fx - j, 0, 1)[:, None], np.clip(fy - i, 0, 1)[:, None]
    v = ((F[p, i, j] * (1 - wx) + F[p, i, j + 1] * wx) * (1 - wy)
         + (F[p, i + 1, j] * (1 - wx) + F[p, i + 1, j + 1] * wx) * wy)
    return v, inside


def _turn(v, phi):
    """d heading / ds (rad/m) for compass heading phi: rays bend towards slower (shallower) water."""
    return (np.sin(phi) * v[:, 2] - np.cos(phi) * v[:, 1]) / v[:, 0]


def trace(x, y, h, sx, sy, reach, deep):
    """Trace the fan of rays from (sx, sy) out to sea for every period. A ray ends once it is `reach`
    metres away and deep[period] metres deep (or MAX_REACH away). Returns launch headings (degrees,
    the direction waves arrive from at the start) and, per [period, ray], the heading and depth where
    the ray ended: NaN where it ran aground, looped or left the grid."""
    dx, dy = x[1] - x[0], y[1] - y[0]
    F = _fields(h, dx, dy)
    launch = np.arange(0, 360, FAN)
    p = np.repeat(np.arange(len(PERIODS)), len(launch))
    phi = np.radians(np.tile(launch, len(PERIODS)))
    px, py = np.full(phi.shape, float(sx)), np.full(phi.shape, float(sy))
    end_dir, end_h = np.full(phi.shape, np.nan), np.full(phi.shape, np.nan)
    live = np.arange(len(phi))
    deep = np.asarray(deep, float)
    for _ in range(int(2 * MAX_REACH / STEP)):
        v, inside = _sample(F, p[live], px[live], py[live], x[0], y[0], dx, dy)
        wet = inside & (v[:, 3] >= DRY)
        far = np.hypot(px[live] - sx, py[live] - sy)
        out = ((far >= reach) & (v[:, 3] >= deep[p[live]])) | (far >= MAX_REACH)
        done = live[out & wet]
        end_dir[done], end_h[done] = np.degrees(phi[done]) % 360, v[out & wet, 3]
        keep = ~out & wet
        live, v = live[keep], v[keep]
        if not len(live):
            break
        a = phi[live]
        mx, my = px[live] + 0.5 * STEP * np.sin(a), py[live] + 0.5 * STEP * np.cos(a)
        ma = a + 0.5 * STEP * _turn(v, a)
        vm, _ = _sample(F, p[live], mx, my, x[0], y[0], dx, dy)
        px[live] += STEP * np.sin(ma)
        py[live] += STEP * np.cos(ma)
        phi[live] = a + STEP * _turn(vm, ma)
    shape = (len(PERIODS), len(launch))
    return launch, end_dir.reshape(shape), end_h.reshape(shape)


def coefficients(launch, end_dir, end_h, depth=TARGET_DEPTH, spread_of=spread):
    """(K, nearshore direction), each [offshore direction, period], from traced rays."""
    K = np.zeros((len(DIRS), len(PERIODS)))
    near = np.full(K.shape, np.nan)
    s, c = np.sin(np.radians(launch)), np.cos(np.radians(launch))
    for p, T in enumerate(PERIODS):
        cn, cgn = _speeds(depth, T)
        co, cgo = _speeds(np.nan_to_num(end_h[p], nan=1.0), T)
        gain = np.where(np.isfinite(end_dir[p]), co * cgo / (cn * cgn), 0.0)
        sd = spread_of(T)
        diff = (np.nan_to_num(end_dir[p])[None, :] - DIRS[:, None] + 180) % 360 - 180
        w = np.exp(-0.5 * (diff / sd) ** 2) / (sd * np.sqrt(2 * np.pi)) * FAN * gain[None, :]
        K[:, p] = np.sqrt(w.sum(1))
        ok = K[:, p] > 0.02
        near[ok, p] = np.degrees(np.arctan2(w[ok] @ s, w[ok] @ c)) % 360
    return K, near


def table(lat, lon, facing, cell, at=None, depth=TARGET_DEPTH, spread_of=spread):
    """The transfer table for one point, as a dict: K and dir [direction, period] and where the rays
    started. at=(lat, lon) starts them nearest that spot instead of the break (validate() uses CDIP's)."""
    x, y, h = _grid(lat, lon)
    px, py = _xy(*at, lat, lon) if at else (0.0, 0.0)
    sx, sy = _start(x, y, h, px, py, facing, depth)
    cx, cy = _xy(cell[0], cell[1], lat, lon)
    reach = max(float(np.hypot(cx - sx, cy - sy)), MIN_REACH)
    # Rays end as deep as the model cell, where its waves apply, or in deep water for their period.
    deep = np.clip(_depth_at(x, y, h, cx, cy), depth, G * PERIODS ** 2 / (4 * np.pi))
    K, near = coefficients(*trace(x, y, h, sx, sy, reach, deep), depth, spread_of)
    window = np.abs((DIRS - facing + 180) % 360 - 180) <= config.EXPOSURE
    return dict(lat=lat, lon=lon, facing=facing, cell=[float(v) for v in cell], depth=depth,
                window_k=round(float(K[window][:, PERIODS == 12].max()), 3),
                start=[round(lat + sy / M_PER_DEG, 5), round(lon + sx / (M_PER_DEG * np.cos(np.radians(lat))), 5)],
                start_offset_m=round(float(np.hypot(sx - px, sy - py))), reach_km=round(reach / 1000, 2),
                dirs=DIRS.tolist(), periods=PERIODS.tolist(), K=np.round(K, 3).tolist(),
                dir=[[None if np.isnan(d) else round(float(d), 1) for d in row] for row in near])


def _cell(lat, lon):
    with xr.open_dataset(points._waves_path(lat, lon)) as w:
        return tuple(float(v) for v in w.attrs["cell"])


def _path(lat, lon):
    return OUT / f"{points._stem(lat, lon)}.json"


def build(workers=4, force=False):
    """Fetch bathymetry and compute the transfer table for every point with a model cell, into
    data/refract/ (kept unless force)."""
    pts = {(la, lo): f for _, _, _, la, lo, f in points.all_points() if points._waves_path(la, lo).exists()}
    t0 = time.time()
    with ThreadPoolExecutor(workers) as pool:
        list(pool.map(lambda p: bathymetry(*p), pts))
    print(f"bathymetry for {len(pts)} points: {time.time() - t0:.0f} s", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    t0, done = time.time(), 0
    for (la, lo), f in pts.items():
        if _path(la, lo).exists() and not force:
            continue
        try:
            _path(la, lo).write_text(json.dumps(table(la, lo, f, _cell(la, lo))))
            done += 1
        except Exception as e:  # keep going; a re-run picks up what failed
            print(f"  {la}, {lo} failed: {e}", flush=True)
    transfer.cache_clear()
    print(f"tables for {done} points: {time.time() - t0:.0f} s", flush=True)


def _apply(t, max_k=None):
    """A function (hs, tp, dir) offshore -> (hs, dir) nearshore from a table dict, with K capped at max_k."""
    K = np.minimum(np.array(t["K"], float), max_k or np.inf)
    near = np.array([[np.nan if d is None else d for d in row] for row in t["dir"]], float)
    E = K * K  # direction vectors weighted by energy interpolate smoothly into a shadow
    u, v = np.nan_to_num(E * np.sin(np.radians(near))), np.nan_to_num(E * np.cos(np.radians(near)))
    dirs, periods = np.array(t["dirs"]), np.array(t["periods"])
    step = dirs[1] - dirs[0]

    def f(hs, tp, dr):
        hs, tp, dr = (np.asarray(a, float) for a in np.broadcast_arrays(hs, tp, dr))
        ok = np.isfinite(hs) & np.isfinite(tp) & np.isfinite(dr)
        fi = np.where(ok, dr % 360, 0) / step
        i0 = np.floor(fi).astype(int) % len(dirs)
        i1, wi = (i0 + 1) % len(dirs), fi - np.floor(fi)
        fj = np.interp(np.where(ok, tp, periods[0]), periods, np.arange(len(periods)))
        j0 = np.minimum(np.floor(fj).astype(int), len(periods) - 2)
        wj = fj - j0
        bil = lambda A: ((A[i0, j0] * (1 - wi) + A[i1, j0] * wi) * (1 - wj)
                         + (A[i0, j0 + 1] * (1 - wi) + A[i1, j0 + 1] * wi) * wj)
        k, e = bil(K), bil(E)
        d = np.where(e > 1e-4, np.degrees(np.arctan2(bil(u), bil(v))) % 360, np.nan)
        return np.where(ok, hs * k, np.nan), np.where(ok, d, np.nan)

    return f


@lru_cache(maxsize=None)
def transfer(lat, lon, max_k=None):
    """(hs, tp, dir) offshore at the point's model cell -> (hs, dir) at TARGET_DEPTH off the break,
    as a vectorised function, or None until build() has made its table. Periods outside PERIODS
    take the nearest column. max_k=1 keeps the losses (shadows, spreading) but not the gains
    (shoaling, focusing), which hang on bathymetric detail ETOPO lacks; see validate(). Also None
    where the table blocks the break's whole window (MIN_WINDOW_K), so callers keep the raw waves."""
    if not _path(lat, lon).exists():
        return None
    t = json.loads(_path(lat, lon).read_text())
    return _apply(t, max_k) if t["window_k"] >= MIN_WINDOW_K else None


# --- Validation against CDIP MOP in California ----------------------------------------------

def _offshore(lat, lon):
    """Copernicus swell trains at a point's mornings, before points.mornings() swaps in CDIP."""
    from .inland import TP_RATIO

    with points._open([points._waves_path(lat, lon, True), points._waves_path(lat, lon)]) as w:
        idx, days = points._local_mornings(w.time.values, lon, 3)
        g = lambda v: w[v].values[idx].astype(float)
        return dict(day=days, hs1=g("VHM0_SW1"), tp1=g("VTM01_SW1") * TP_RATIO, dir1=g("VMDR_SW1"),
                    hs2=g("VHM0_SW2"), tp2=g("VTM01_SW2") * TP_RATIO, dir2=g("VMDR_SW2"))


def _stats(s, facing):
    """For one source's two trains: swell Hs (combined), dominant direction, the share of mornings
    with a train of 4 ft over 9 s (big; big_win also within the break's window) and the mornings
    whose swell passes the good+ test (rating.swell_ok)."""
    from . import rating

    h1, h2 = np.nan_to_num(s["hs1"]), np.nan_to_num(s["hs2"])
    dom = np.where(h2 > h1, s["dir2"], s["dir1"])
    big = lambda k, win: ((s[f"hs{k}"] >= config.MIN_SWELL_HEIGHT) & (s[f"tp{k}"] > config.MIN_SWELL_PERIOD)
                          & (rating.exposed(s[f"dir{k}"], facing) if win else True))
    ok = rating.swell_ok(s["hs1"], s["tp1"], s["dir1"], facing) | rating.swell_ok(s["hs2"], s["tp2"], s["dir2"], facing)
    return dict(hs=np.hypot(h1, h2), dir=dom, big=100 * (big(1, False) | big(2, False)).mean(),
                big_win=100 * (big(1, True) | big(2, True)).mean(), ok=ok)


def validate(max_k=None, at_site=False, spread_of=spread):
    """Compare raw offshore and refracted Copernicus swell with CDIP MOP at each California break,
    over the same mornings; prints a table and returns its rows. "agree" is the share of mornings on
    which the good+ swell test gives CDIP's answer. at_site=True starts the rays at the CDIP site
    and depth (on the fly) rather than using the break's saved table."""
    from . import cdip

    rows = []
    for kind, name, region, lat, lon, facing in points.all_points():
        near = cdip.mornings(lat, lon)
        if near is None:
            continue
        if at_site or spread_of is not spread:
            with xr.open_dataset(cdip._path(near["site"])) as d:
                site = (d.attrs["lat"], d.attrs["lon"]) if at_site else None
                depth = d.attrs["depth"] if at_site else TARGET_DEPTH
            f = _apply(table(lat, lon, facing, _cell(lat, lon), site, depth, spread_of), max_k)
        else:
            f = transfer(lat, lon, max_k)
        off = _offshore(lat, lon)
        _, io, ic = np.intersect1d(off["day"], near["day"], return_indices=True)
        o = {k: v[io] for k, v in off.items() if k != "day"}
        c = {k: near[k][ic] for k in ("hs1", "tp1", "dir1", "hs2", "tp2", "dir2")}
        keep = np.isfinite(o["hs1"]) & np.isfinite(c["hs1"])
        o, c = ({k: v[keep] for k, v in d.items()} for d in (o, c))
        r = dict(o)
        r["hs1"], r["dir1"] = f(o["hs1"], o["tp1"], o["dir1"])
        r["hs2"], r["dir2"] = f(o["hs2"], o["tp2"], o["dir2"])
        sc, so, sr = (_stats(s, facing) for s in (c, o, r))
        swell = c["tp1"] > config.MIN_SWELL_PERIOD
        row = dict(name=name, site=near["site"], n=int(keep.sum()), cdip=sc["big"], cdip_win=sc["big_win"])
        for tag, s in (("raw", so), ("ref", sr)):
            err = s["hs"] - sc["hs"]
            ddir = np.abs((s["dir"] - sc["dir"] + 180) % 360 - 180)[swell]
            row.update({f"{tag}_bias": err.mean(), f"{tag}_rmse": np.sqrt((err ** 2).mean()),
                        f"{tag}_r": np.corrcoef(s["hs"], sc["hs"])[0, 1], f"{tag}_ddir": np.nanmedian(ddir),
                        tag: s["big"], f"{tag}_win": s["big_win"], f"{tag}_agree": 100 * (s["ok"] == sc["ok"]).mean()})
        rows.append(row)
    head = ("break", "site", "n", "bias raw", "bias ref", "rmse raw", "rmse ref", "r raw", "r ref",
            "ddir raw", "ddir ref", "%4ft CDIP", "raw", "raw(win)", "ref", "agree raw", "agree ref")
    print(("{:18} {:6} {:>5}" + " {:>8}" * 14).format(*head))
    for w in rows:
        print(("{:18} {:6} {:5d}" + " {:8.2f}" * 6 + " {:8.0f}" * 2 + " {:8.1f}" * 6).format(
            w["name"], w["site"], w["n"], w["raw_bias"], w["ref_bias"], w["raw_rmse"], w["ref_rmse"],
            w["raw_r"], w["ref_r"], w["raw_ddir"], w["ref_ddir"], w["cdip"], w["raw"], w["raw_win"], w["ref"],
            w["raw_agree"], w["ref_agree"]))
    return rows
