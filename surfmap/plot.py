"""Render one PDF page per month: world contour map, enclosed-sea insets, best-breaks table."""

import calendar

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shp
import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import BoundaryNorm, ListedColormap
from shapely.ops import unary_union

from . import config, inland, rating
from .analysis import LAT, LON, Period

# Red (rarely good) through pale yellow (~the open-ocean median) to green (usually good).
# The green end is much darker than the red end, so red-green colour-blind readers can
# still tell the extremes apart by lightness.
LEVELS = [0, 2, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100]
COLORS = ["#d73027", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", "#d9ef8b",
          "#a6d96a", "#66bd63", "#1a9850", "#006837", "#00441b"]
CMAP, NORM = ListedColormap(COLORS), BoundaryNorm(LEVELS, len(COLORS))

INK, INK_2, INK_3 = "#1a1a19", "#4a4a47", "#7a7a76"
LAND, COAST, NODATA, HATCH = "#e4e3de", "#9a9994", "#d0cfca", "#a9a8a2"
PC = ccrs.PlateCarree()
MAP_CENTER = 10  # central longitude of the world maps
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]

INSETS = [  # title, basins it shows, extent (lon0, lon1, lat0, lat1), regional pages whose best coasts it labels
    ("Great Lakes", ("great_lakes",), (-94.2, -75.5, 40.8, 49.3), (("Lake Michigan", 2), ("Lake Superior", 2))),
    ("Mediterranean & Black Sea", ("mediterranean", "black_sea"), (-6.5, 42.5, 30.0, 47.5),
     (("Mediterranean", 4), ("Greece", 1), ("Black Sea", 2))),
    ("Caspian Sea", ("caspian",), (45.8, 55.5, 36.2, 47.5), (("Caspian Sea", 3),)),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
    "hatch.color": HATCH, "hatch.linewidth": 0.5, "pdf.fonttype": 42,
})


def _land():
    """Natural Earth land with the Great Lakes cut out (the Caspian is already a hole)."""
    land = unary_union(list(shp.Reader(shp.natural_earth("50m", "physical", "land")).geometries()))
    return cfeature.ShapelyFeature([land.difference(inland.geometries()["great_lakes"])], PC)


def _draw(ax, field, nodata, land, rows=slice(None), cols=slice(None), center=None):
    """Contour the field onto ax. The grid is projected first and contoured in map coordinates:
    reprojecting finished contours breaks on huge globe-spanning bands (January's 0-2% band came
    out in the wrong colour). For a global map, pass its central longitude so the grid starts at
    the map's edge and no contour straddles the seam."""
    lo, la = LON[cols], LAT[rows]
    f, nd = field[rows][:, cols], nodata[rows][:, cols].astype(float)
    if center is not None:
        k = int(np.searchsorted(lo, center - 180))
        lo = np.concatenate([lo[k:], lo[:k] + 360])
        f, nd = (np.concatenate([a[:, k:], a[:, :k]], axis=1) for a in (f, nd))
    x, y = np.meshgrid(lo, la)
    cs = ax.contourf(x, y, f, levels=LEVELS, cmap=CMAP, norm=NORM, transform=PC, transform_first=True, zorder=1)
    cs.set_rasterized(True)
    hatch = ax.contourf(x, y, nd, levels=[0.5, 1.5], colors=[NODATA], hatches=["////"], transform=PC,
                        transform_first=True, zorder=1.5)
    hatch.set_rasterized(True)
    # Rasterized: as vectors the 1:50m outlines make each page ~4 MB instead of ~1 MB.
    ax.add_feature(land, facecolor=LAND, edgecolor="none", zorder=2, rasterized=True)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=COAST, linewidth=0.3, zorder=3,
                   rasterized=True)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="none", edgecolor=COAST,
                   linewidth=0.3, zorder=3, rasterized=True)
    ax.spines["geo"].set_edgecolor(COAST)
    ax.spines["geo"].set_linewidth(0.5)


def _swell_arrows(ax, lon, lat, u, v):
    q = ax.quiver(lon, lat, u, v, transform=PC, pivot="middle", color=INK, alpha=0.8,
                  scale=55, width=0.0013, headwidth=4.5, headlength=5, headaxislength=4.5, zorder=4)
    q.set_path_effects([pe.withStroke(linewidth=1.2, foreground="white", alpha=0.6)])
    return q


def _label_spots(ax, spots, n=16, min_sep=7.0):
    """Dots for every break; names for the best few, skipping ones crowded by a better neighbour."""
    lons = [s["lon"] for s in spots]
    lats = [s["lat"] for s in spots]
    ax.scatter(lons, lats, s=9, c=INK, edgecolors="white", linewidths=0.6, transform=PC, zorder=5)
    placed = []
    for rank, s in enumerate(spots[:40], 1):
        if len(placed) >= n:
            break
        if any(np.hypot(s["lon"] - a, s["lat"] - b) < min_sep for a, b in placed):
            continue
        placed.append((s["lon"], s["lat"]))
        ax.text(s["lon"] + 1.2, s["lat"] + 0.6, f"{s['name']} {s['pct']:.0f}%", fontsize=6.8,
                color=INK, transform=PC, zorder=6, path_effects=HALO)


def _stretch_markers(ax, rows):
    """Numbered diamonds for the under-the-radar stretches; numbers follow stretches.py."""
    ax.scatter([r["lon"] for r in rows], [r["lat"] for r in rows], s=22, marker="D", c="white",
               edgecolors=INK, linewidths=0.8, transform=PC, zorder=6)
    for r in rows:
        ax.text(r["lon"] - 1.2, r["lat"], str(r["num"]), fontsize=6.3, weight="bold", ha="right",
                va="center", color=INK, transform=PC, zorder=7, path_effects=HALO)


TABLE_SEP_KM = 161  # 100 miles: one entry per stretch of coast, the best-ranked one


def _km(a, b):
    """Great-circle distance between two spots."""
    la1, lo1, la2, lo2 = map(np.radians, (a["lat"], a["lon"], b["lat"], b["lon"]))
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371 * np.arcsin(np.sqrt(h))


def _distinct(spots, n):
    """The best n spots, skipping any within TABLE_SEP_KM of one already chosen."""
    chosen = []
    for s in spots:
        if all(_km(s, c) >= TABLE_SEP_KM for c in chosen):
            chosen.append(s)
            if len(chosen) == n:
                break
    return chosen


def _table(fig, rect, spots, stretches, title, breakout_page):
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.text(0, 1, title, fontsize=12, weight="bold", va="top")
    ax.text(0, 0.968, "% of mornings good+ at the break: swell inside its window, wind\n"
            "light or offshore; best per 100 miles; then typical morning wind",
            fontsize=8, color=INK_2, va="top", linespacing=1.3)
    top = _distinct([s for s in spots if s["basin"] is None], 15)
    y, step = 0.905, 0.0375
    for i, s in enumerate(top, 1):
        ax.text(0.0, y, f"{i}", fontsize=8.5, color=INK_3, va="top")
        ax.text(0.07, y, s["name"], fontsize=8.5, weight="bold", va="top")
        ax.text(0.07, y - 0.017, " · ".join(filter(None, (s["region"], _wind_text(s["wind"])))),
                fontsize=7, color=INK_2, va="top")
        ax.text(1.0, y, f"{s['pct']:.0f}%", fontsize=9, weight="bold", va="top", ha="right")
        y -= step

    y -= 0.012
    ax.text(0, y, "Under-the-radar stretches", fontsize=11, weight="bold", va="top")
    ax.text(0, y - 0.03, f"◆ numbered on the map; month by month on page {breakout_page}",
            fontsize=8, color=INK_2, va="top")
    y -= 0.07
    for s in sorted(stretches, key=lambda r: -r["pct"])[:7]:
        ax.text(0.0, y, f"◆{s['num']}", fontsize=8.5, color=INK_3, va="top")
        ax.text(0.09, y, s["name"], fontsize=8.5, weight="bold", va="top")
        ax.text(0.09, y - 0.017, " · ".join(filter(None, (s["region"], _wind_text(s["wind"])))),
                fontsize=7, color=INK_2, va="top")
        ax.text(1.0, y, f"{s['pct']:.0f}%", fontsize=9, weight="bold", va="top", ha="right")
        y -= step


def _inset_labels(sources, m):
    """The best coasts from each regional page (fixed across pages), with this period's % under the
    relaxed enclosed-sea standard the inset is contoured with."""
    from . import inland, regions

    ratio = m.tp_ratio or inland.TP_RATIO
    labels, seen = [], set()
    for region, n in sources:
        values = regions.coast_pct(region, m.months, regions.LAKE, ratio)
        for name in regions.best(region, n):
            if name in values and name not in seen:
                lat, lon, pct = values[name]
                labels.append(dict(name=name, lat=lat, lon=lon, pct=pct))
                seen.add(name)
    return labels


def _inset(fig, rect, title, basins, extent, field, nodata, land, labels, ready):
    lon0, lon1, lat0, lat1 = extent
    ax = fig.add_axes(rect, projection=ccrs.Mercator(central_longitude=(lon0 + lon1) / 2))
    ax.set_extent(extent, crs=PC)
    pad = 3
    rows = (LAT > lat0 - pad) & (LAT < lat1 + pad)
    cols = (LON > lon0 - pad) & (LON < lon1 + pad)
    _draw(ax, field, nodata, land, rows, cols)
    ax.scatter([s["lon"] for s in labels], [s["lat"] for s in labels], s=10, c=INK,
               edgecolors="white", linewidths=0.6, transform=PC, zorder=5)
    for s in labels:
        # Close neighbours split sideways (west one to the left, east one to the right); others sit above.
        near = [o for o in labels if o is not s and abs(o["lon"] - s["lon"]) < 0.25 * (lon1 - lon0)
                and abs(o["lat"] - s["lat"]) < 0.12 * (lat1 - lat0)]
        if near:
            west = all(s["lon"] < o["lon"] for o in near)
            x, y, ha, va = s["lon"] + (-1 if west else 1) * 0.015 * (lon1 - lon0), s["lat"], "right" if west else "left", "center"
        else:
            x, y, ha, va = s["lon"], s["lat"] + 0.036 * (lat1 - lat0), "center", "bottom"
        ax.text(x, y, f"{s['name']} {s['pct']:.0f}%", fontsize=6.5, ha=ha, va=va, transform=PC, zorder=6,
                path_effects=HALO)
    sep = "\n" if rect[2] < 0.15 else "  \u00b7  "  # narrow insets put the standard on its own line
    ax.set_title(f"{title}{sep}good = 4 ft @ 6 s", fontsize=9.5, weight="bold", loc="left", pad=4)
    if any(b not in ready for b in basins):
        ax.text(0.5, 0.5, "Open-Meteo data is\nstill downloading", transform=ax.transAxes,
                ha="center", va="center", fontsize=8.5, color=INK_2, zorder=7,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor=COAST, linewidth=0.5))


def _inset_rects(fig, bottom, height, left=0.02, right=0.98):
    """Inset boxes sized to each map's projected aspect, spread with equal gaps across the page."""
    fw, fh = fig.get_size_inches()
    widths = []
    for _, _, (lon0, lon1, lat0, lat1), _ in INSETS:
        proj = ccrs.Mercator(central_longitude=(lon0 + lon1) / 2)
        (x0, y0, _), (x1, y1, _) = proj.transform_points(PC, np.array([lon0, lon1]), np.array([lat0, lat1]))
        widths.append(height * fh * (x1 - x0) / (y1 - y0) / fw)
    gap = (right - left - sum(widths)) / (len(widths) - 1)
    x, rects = left, []
    for w in widths:
        rects.append([x, bottom, w, height])
        x += w + gap
    return rects


def _colorbar(fig, rect, mappable):
    cax = fig.add_axes(rect)
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", spacing="uniform", ticks=LEVELS,
                      format=matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}%"))
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, length=0, colors=INK_2)
    cax.set_title("Share of mornings with good-or-better surf", fontsize=8.5, color=INK_2, loc="left", pad=3)
    fig.patches.append(matplotlib.patches.Rectangle(
        (rect[0] + rect[2] + 0.02, rect[1]), 0.022, rect[3], transform=fig.transFigure,
        facecolor=NODATA, hatch="////", edgecolor=HATCH, linewidth=0))
    fig.text(rect[0] + rect[2] + 0.046, rect[1] + rect[3] / 2, "No model data (sea ice, or a\nsea neither source covers)",
             fontsize=7.5, color=INK_2, va="center", linespacing=1.2)


def page(pdf, m, land, breakout_page):
    """Render one Period as a PDF page."""
    field, nodata = m.field()
    spots = m.spots()
    stretches = m.stretches()

    fig = plt.figure(figsize=(17, 11))
    fig.patch.set_facecolor("white")
    fig.text(0.02, 0.965, f"{m.title}: where the surf is good", fontsize=24, weight="bold", va="top")
    fig.text(0.02, 0.925, f"Share of mornings with good-or-better surf in the wave-model record, "
             f"{m.describe()}", fontsize=11.5, color=INK_2, va="top")

    ax = fig.add_axes([0.015, 0.375, 0.73, 0.54], projection=ccrs.Robinson(central_longitude=MAP_CENTER))
    ax.set_global()
    _draw(ax, field, nodata, land, center=MAP_CENTER)
    arrows = _swell_arrows(ax, *m.swell_arrows())
    _label_spots(ax, [s for s in spots if s["basin"] is None])
    _stretch_markers(ax, stretches)
    _colorbar(fig, [0.06, 0.345, 0.46, 0.013], ax.collections[0])
    ax.quiverkey(arrows, 0.715, 0.3515, 1, "Predominant swell direction (power-weighted;\n"
                 "omitted where swell comes from mixed directions)", coordinates="figure",
                 labelpos="E", labelsep=0.08, fontproperties=dict(size=7.5), labelcolor=INK_2)

    _table(fig, [0.77, 0.36, 0.215, 0.57], spots, stretches, m.table_title, breakout_page)

    for (title, basins, extent, sources), rect in zip(INSETS, _inset_rects(fig, bottom=0.068, height=0.235)):
        _inset(fig, rect, title, basins, extent, field, nodata, land, _inset_labels(sources, m), m.ready)

    fig.text(0.02, 0.045, rating.rule_text(), fontsize=7.5, color=INK_2, va="top")
    fig.text(0.02, 0.028,
             "Ocean: NOAA WAVEWATCH III global model (Univ. of Hawaiʻi run, PacIOOS ERDDAP), 1° sampling. "
             "Enclosed seas (MFWAM via Open-Meteo and Copernicus Marine; 0.75° lakes, 1.5–2° seas) use a relaxed "
             "standard: total waves 4–10 ft at a dominant period of at least 6 s.",
             fontsize=7.5, color=INK_2, va="top")
    fig.text(0.02, 0.011,
             "Map colours: offshore sea state only (no bathymetry, sheltering or local wind). Marked breaks and "
             "\u25c6 stretches: either swell train at each break's ~8 km MFWAM cell (CDIP nearshore in California), "
             f"with ERA5 morning wind \u2264 {config.MAX_WIND_KT:g} kt or offshore instead of the wind-sea test.",
             fontsize=7.5, color=INK_3, va="top")
    pdf.savefig(fig, dpi=220)
    return fig


def _compass(deg, degrees=True):
    if np.isnan(deg):
        return "–"
    points = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    name = points[int((deg + 11.25) // 22.5) % 16]
    return f"{name} {deg:.0f}°" if degrees else name


def _wind_text(w):
    """Typical morning wind as '7 kt SSE', or '' where there is none."""
    return f"{w[0]:.0f} kt {_compass(w[1], degrees=False)}" if w else ""


def _ink_on(color):
    """Dark or white text, whichever reads on this fill."""
    r, g, b = (c if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in matplotlib.colors.to_rgb(color))
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.25 else INK


def _peak(vals):
    """Best three consecutive calendar months: (label such as 'Jun\u2013Aug', mean %)."""
    v = np.asarray(vals, float)
    means = [np.nanmean(np.roll(v, -s)[:3]) for s in range(12)]
    s = int(np.nanargmax(means))
    return f"{calendar.month_abbr[s + 1]}\u2013{calendar.month_abbr[(s + 2) % 12 + 1]}", means[s]


TIDE_NOTE = "Tide range = typical to big-tide daily range (EOT20, DGFI-TUM)."


def _tide_text(name):
    """'5\u20137 ft' from tides.py; 'none' where the model has no ocean tide (Great Lakes, Caspian)."""
    from . import tides

    return tides.text(name, rating.FT) or "none"


WIND_NOTE = ("Morning wind = median speed and most common direction over all mornings of the period, "
             "ERA5 10 m wind just offshore.")


def _breakout_page(pdf, title, subtitle, kind, rows, notes):
    """Month-by-month table. Each row needs tag, name, region, months (12 % for Jan..Dec), pct, the
    median good morning (hs, hb, tp, dir) and wind."""
    fig = plt.figure(figsize=(17, 11))
    fig.patch.set_facecolor("white")
    fig.text(0.02, 0.965, title, fontsize=24, weight="bold", va="top")
    fig.text(0.02, 0.925, subtitle, fontsize=11.5, color=INK_2, va="top", linespacing=1.4)

    ax = fig.add_axes([0.02, 0.075, 0.96, 0.8])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    h = 1 / (len(rows) + 1.5)
    x_month, w_month = 0.30, 0.028
    x_year, x_worst, x_best = 0.655, 0.697, 0.72
    x_tide, x_swell, x_wind = 0.78, 0.828, 0.94
    y = 1 - h / 2
    for x, label, ha in [(0.0, "#", "left"), (0.025, kind, "left"), (0.165, "Region", "left"),
                         (x_year, "Year", "center"), (x_worst, "Worst", "center"), (x_best, "Best 3 months", "left"),
                         (x_tide, "Tide range", "left"), (x_swell, "Typical good morning", "left"),
                         (x_wind, "Morning wind", "left")]:
        ax.text(x, y, label, fontsize=8.5, weight="bold", color=INK_2, ha=ha, va="center")
    for k, mon in enumerate("JFMAMJJASOND"):
        ax.text(x_month + (k + 0.5) * w_month, y, mon, fontsize=8.5, weight="bold", color=INK_2,
                ha="center", va="center")

    size = 8 if len(rows) > 36 else 8.5
    for r in rows:
        y -= h
        ax.text(0.0, y, r["tag"], fontsize=size, color=INK_3, va="center")
        ax.text(0.025, y, r["name"], fontsize=size, weight="bold", va="center")
        ax.text(0.165, y, r["region"], fontsize=size - 0.5, color=INK_2, va="center")
        for k, v in enumerate(r["months"]):
            color = CMAP(NORM(v)) if np.isfinite(v) else NODATA
            ax.add_patch(matplotlib.patches.Rectangle((x_month + k * w_month + 0.001, y - h * 0.42),
                                                      w_month - 0.002, h * 0.84, facecolor=color, linewidth=0))
            ax.text(x_month + (k + 0.5) * w_month, y, f"{v:.0f}" if np.isfinite(v) else "\u2013",
                    fontsize=size - 1, ha="center", va="center", color=_ink_on(color))
        best, best_pct = _peak(r["months"])
        ax.text(x_year, y, f"{r['pct']:.0f}%", fontsize=size + 0.5, weight="bold", ha="center", va="center")
        ax.text(x_worst, y, f"{np.nanmin(r['months']):.0f}%", fontsize=size, ha="center", va="center")
        ax.text(x_best, y, f"{best} {best_pct:.0f}%", fontsize=size, va="center")
        good = (f"{r['hs'] * rating.FT:.0f} ft @ {r['tp']:.0f} s {_compass(r['dir'], degrees=False)}, "
                f"surf ~{r['hb'] * rating.FT:.0f} ft") if np.isfinite(r["hs"]) else "\u2013"
        ax.text(x_tide, y, _tide_text(r["name"]), fontsize=size, color=INK_2, va="center")
        ax.text(x_swell, y, good, fontsize=size, color=INK_2, va="center")
        ax.text(x_wind, y, _wind_text(r["wind"]) or "\u2013", fontsize=size, color=INK_2, va="center")

    fig.text(0.02, 0.06, "Month cells use the map colours (% of mornings good+). Typical good morning = medians "
             "over the good+ mornings: swell height @ peak period and the direction it comes from, and estimated "
             "breaking surf (Komar\u2013Gaughan). " + TIDE_NOTE + "\n" + WIND_NOTE, fontsize=7.5, color=INK_2, va="top",
             linespacing=1.5)
    fig.text(0.02, 0.026, notes, fontsize=7.5, color=INK_3, va="top")
    pdf.savefig(fig, dpi=220)
    return fig


def _monthly(kind, key):
    """{key value: [% good+ for Jan..Dec]} for 'spots' or 'stretches'."""
    out = {}
    for m in range(1, 13):
        for r in getattr(Period((m,), "", ""), kind)():
            out.setdefault(r[key], [np.nan] * 12)[m - 1] = r["pct"]
    return out


def stretch_page(pdf):
    """One row per discovered stretch, in stretches.py order."""
    from . import discover

    monthly = _monthly("stretches", "num")
    rows = sorted((dict(r, tag=f"\u25c6{r['num']}", months=monthly[r["num"]])
                   for r in Period(tuple(range(1, 13)), "", "").stretches()), key=lambda r: -r["pct"])
    subtitle = (f"Coast where mornings are often good+ with swell aimed at the shore (\u2265{discover.MIN_YEAR}% "
                f"year-round, or \u2265{discover.MIN_PEAK}% over its best three months),\nat least {discover.FAR_KM} "
                f"km from the well-known breaks and within {discover.TOWN_KM} km of a town of "
                f"{discover.TOWN_POP:,}+. Listed by year-round consistency; numbers match the \u25c6 markers on the maps.")
    notes = (f"Stretches were found on the WW3 grid; each is rated at its coast from M\u00e9t\u00e9o-France MFWAM "
             f"(Copernicus Marine, ~8 km), either swell train, with ERA5 morning wind \u2264 {config.MAX_WIND_KT:g} kt "
             f"or offshore. Check local reefs and access before booking a trip.")
    return _breakout_page(pdf, "Under-the-radar stretches, month by month", subtitle, "Stretch", rows, notes)


POPULAR_ROWS = 40


def spot_page(pdf):
    """The named ocean breaks with the strongest three-month season, one per 100 miles."""
    monthly = _monthly("spots", "name")
    rows = [dict(r, months=monthly[r["name"]]) for r in Period(tuple(range(1, 13)), "", "").spots()
            if r["basin"] is None and r["name"] in monthly]
    rows = _distinct(sorted(rows, key=lambda r: -_peak(r["months"])[1]), POPULAR_ROWS)
    rows.sort(key=lambda r: -r["pct"])  # listed by year-round consistency
    for i, r in enumerate(rows, 1):
        r["tag"] = str(i)
    subtitle = (f"The {len(rows)} well-known breaks with the strongest three-month season (best one per 100 miles "
                f"of coast), listed by year-round consistency.\nBreaks in the enclosed seas are on the map insets.")
    notes = (f"Each break is read from M\u00e9t\u00e9o-France MFWAM (Copernicus Marine) at its ~8 km cell: either swell "
             f"train inside its exposure window, with ERA5 morning wind \u2264 {config.MAX_WIND_KT:g} kt or offshore. "
             f"Reefs, points and sheltering inside that cell are not modelled.")
    return _breakout_page(pdf, "Popular breaks, month by month", subtitle, "Break", rows, notes)


def _label_spot(r, rows, extent):
    """Where to put a coast's name: seaward first, then the other sides, keeping clear of the
    other markers and inside the map. Sizes are rough, in degrees of the map extent."""
    lon0, lon1, lat0, lat1 = extent
    gap_x, gap_y = 0.03 * (lon1 - lon0), 0.03 * (lat1 - lat0)
    w, h = 0.0125 * (lon1 - lon0) * len(r["name"]), 0.035 * (lat1 - lat0)
    best = None
    for turn in (0, 180, 90, -90):
        b = np.radians(r["facing"] + turn)
        dx, dy = np.sin(b), np.cos(b)
        ha = "left" if dx > 0.35 else "right" if dx < -0.35 else "center"
        va = "bottom" if dy > 0.35 else "top" if dy < -0.35 else "center"
        x, y = r["lon"] + dx * gap_x, r["lat"] + dy * gap_y
        x0 = {"left": x, "right": x - w, "center": x - w / 2}[ha]
        y0 = {"bottom": y, "top": y - h, "center": y - h / 2}[va]
        inside = lon0 <= x0 and x0 + w <= lon1 and lat0 <= y0 and y0 + h <= lat1
        clear = all(not (x0 - gap_x / 2 <= o["lon"] <= x0 + w + gap_x / 2 and y0 - gap_y / 2 <= o["lat"] <= y0 + h + gap_y / 2)
                    for o in rows if o is not r)
        best = best or (x, y, ha, va)
        if inside and clear:
            return x, y, ha, va
    return best


SPLIT_ROWS = 16  # regions with more coasts than this get a map page and a separate table page


def _new_page(title, subtitle):
    fig = plt.figure(figsize=(17, 11))
    fig.patch.set_facecolor("white")
    fig.text(0.02, 0.965, title, fontsize=24, weight="bold", va="top")
    fig.text(0.02, 0.925, subtitle, fontsize=11.5, color=INK_2, va="top", linespacing=1.4)
    return fig


def _region_subtitle(std, tp_ratio):
    if std.get("ocean"):  # the good+ rule in brief; the map pages' footer has it in full
        return (f"Good+ = swell \u2265 {config.MIN_SWELL_HEIGHT * rating.FT:.0f} ft over {config.MIN_SWELL_PERIOD:g} s, or "
                f"\u2265 {config.LONG_SWELL_HEIGHT * rating.FT:.0f} ft over {config.LONG_SWELL_PERIOD:g} s; breaking under "
                f"{config.MAX_BREAKER * rating.FT:.0f} ft; either swell train, aimed at the shore; morning wind \u2264 "
                f"{config.MAX_WIND_KT:g} kt or offshore; on the ~{config.LOCAL_HOUR:02d}:00 sample.\nRead from "
                f"M\u00e9t\u00e9o-France MFWAM (Copernicus Marine) at each coast's ~8 km cell, with ERA5 wind just offshore.")
    period = (f"dominant period \u2265 {std['min_period']:g} s (taken as {tp_ratio:.2f}\u00d7 the mean period)"
              if std["peak"] else f"mean period \u2265 {std['min_period']:g} s")
    return (f"{std['name']} = total wave height {std['min_hs'] * rating.FT:.0f}\u2013{std['max_hs'] * rating.FT:.0f} ft "
            f"({std['min_hs']:g}\u2013{std['max_hs']:g} m) with {period}, arriving within {config.EXPOSURE}\u00b0 of "
            f"the shore's facing, on the ~{config.LOCAL_HOUR:02d}:00 sample.\nEnclosed seas rarely see groundswell, "
            f"so this regional standard sits alongside ocean good+, the open-ocean standard used on the other pages.")


def _region_map(fig, rect, spec, rows, land):
    """Numbered coasts, coloured by year-round surfable %."""
    lon0, lon1, lat0, lat1 = spec["extent"]
    ax = fig.add_axes(rect, projection=ccrs.Mercator(central_longitude=(lon0 + lon1) / 2))
    ax.set_extent(spec["extent"], crs=PC)
    ax.set_facecolor("#f4f7fa")
    ax.add_feature(land, facecolor=LAND, edgecolor="none", zorder=2, rasterized=True)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=COAST, linewidth=0.4, zorder=3)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="none", edgecolor=COAST, linewidth=0.4, zorder=3)
    ax.spines["geo"].set_edgecolor(COAST)
    for i, r in enumerate(rows, 1):
        color = CMAP(NORM(r["pct"])) if r["complete"] else NODATA
        ax.scatter([r["lon"]], [r["lat"]], s=170, c=[color], edgecolors=INK, linewidths=0.8, transform=PC, zorder=5)
        ax.text(r["lon"], r["lat"], str(i), fontsize=7.5, weight="bold", ha="center", va="center",
                color=_ink_on(color), transform=PC, zorder=6)
        if r["label"]:  # placed by hand, with a leader line back to the marker
            geo = PC._as_mpl_transform(ax)
            ax.annotate(r["name"], xy=(r["lon"], r["lat"]), xytext=r["label"], xycoords=geo, textcoords=geo,
                        fontsize=8, ha="right" if r["label"][0] < r["lon"] - 0.05 else
                        "left" if r["label"][0] > r["lon"] + 0.05 else "center", va="center", zorder=6,
                        path_effects=HALO, arrowprops=dict(arrowstyle="-", color=INK_3, lw=0.6, shrinkB=6))
        else:
            x, y, ha, va = _label_spot(r, rows, spec["extent"])
            ax.text(x, y, r["name"], fontsize=8, ha=ha, va=va, transform=PC, zorder=6, path_effects=HALO)


def _region_season(fig, rect, rows, std):
    """The region's average surfable % by calendar month."""
    months = np.array([r["months"] for r in rows])
    avg = [np.mean(months[:, k]) if np.isfinite(months[:, k]).all() else np.nan for k in range(12)]
    bx = fig.add_axes(rect)
    for k, v in enumerate(avg):
        if np.isfinite(v):
            bx.bar(k, v, width=0.72, color=CMAP(NORM(v)), zorder=2)
            bx.text(k, v + 0.6, f"{v:.0f}%", ha="center", va="bottom", fontsize=8.5, color=INK)
        else:
            bx.text(k, 0.5, "pending", ha="center", va="bottom", fontsize=7, color=INK_3, rotation=90)
    bx.set_xticks(range(12), [calendar.month_abbr[m] for m in range(1, 13)])
    bx.set_title(f"Average of these coasts: % of mornings with {std['name'].lower()}", fontsize=10.5,
                 weight="bold", loc="left", pad=8)
    bx.spines[["top", "right", "left"]].set_visible(False)
    bx.spines["bottom"].set_color(COAST)
    bx.tick_params(axis="y", left=False, labelleft=False)
    bx.tick_params(axis="x", length=0, labelsize=9, colors=INK_2)
    bx.set_ylim(0, max(10, np.nanmax(avg) * 1.18) if np.isfinite(avg).any() else 10)


def _region_table(fig, rect, rows, std):
    tx = fig.add_axes(rect)
    tx.axis("off")
    tx.set_xlim(0, 1)
    tx.set_ylim(0, 1)
    h = 1 / (len(rows) + 1.3)
    x_month, w_month = 0.25, 0.03
    x_year, x_best, x_ocean, x_tide, x_good, x_wind = 0.635, 0.665, 0.765, 0.795, 0.84, 0.945
    y = 1 - h / 2
    for x, label, ha in [(0.0, "#", "left"), (0.022, "Coast", "left"), (0.12, "Area", "left"),
                         (x_year, "Year", "center"), (x_best, "Best 3 months", "left"),
                         (x_ocean, "Ocean good+", "center"), (x_tide, "Tide range", "left"),
                         (x_good, f"Typical {std['noun']} morning", "left"),
                         (x_wind, "Morning wind", "left")]:
        tx.text(x, y, label, fontsize=8.5, weight="bold", color=INK_2, ha=ha, va="center")
    for k, mon in enumerate("JFMAMJJASOND"):
        tx.text(x_month + (k + 0.5) * w_month, y, mon, fontsize=8.5, weight="bold", color=INK_2,
                ha="center", va="center")
    for i, r in enumerate(rows, 1):
        y -= h
        tx.text(0.0, y, str(i), fontsize=8.5, color=INK_3, va="center")
        tx.text(0.022, y, r["name"], fontsize=8.5, weight="bold", va="center")
        tx.text(0.12, y, r["region"], fontsize=8, color=INK_2, va="center")
        for k, v in enumerate(r["months"]):
            color = CMAP(NORM(v)) if np.isfinite(v) else NODATA
            tx.add_patch(matplotlib.patches.Rectangle((x_month + k * w_month + 0.001, y - h * 0.42),
                                                      w_month - 0.002, h * 0.84, facecolor=color, linewidth=0))
            tx.text(x_month + (k + 0.5) * w_month, y, f"{v:.0f}" if np.isfinite(v) else "\u2013",
                    fontsize=7.5, ha="center", va="center", color=_ink_on(color))
        if r["complete"]:
            best, best_pct = _peak(r["months"])
            tx.text(x_year, y, f"{r['pct']:.0f}%", fontsize=9, weight="bold", ha="center", va="center")
            tx.text(x_best, y, f"{best} {best_pct:.0f}%", fontsize=8.5, va="center")
            tx.text(x_ocean, y, f"{r['ocean']:.0f}%", fontsize=8.5, ha="center", va="center")
        else:
            tx.text(x_year, y, "pending", fontsize=8, color=INK_3, ha="center", va="center")
        good = (f"{r['hs'] * rating.FT:.0f} ft @ {r['tp']:.0f} s {_compass(r['dir'], degrees=False)}"
                if np.isfinite(r["hs"]) else "\u2013")
        tx.text(x_tide, y, _tide_text(r["name"]), fontsize=8.5, color=INK_2, va="center")
        tx.text(x_good, y, good, fontsize=8.5, color=INK_2, va="center")
        tx.text(x_wind, y, _wind_text(r["wind"]) or "\u2013", fontsize=8.5, color=INK_2, va="center")


def _region_footer(fig, rows, std):
    pending = "" if all(r["complete"] for r in rows) else \
        " Months marked \u2013 are still downloading."
    fig.text(0.02, 0.05, "M\u00e9t\u00e9o-France MFWAM (Copernicus Marine) at each coast's ~8 km cell, a few km "
             f"offshore; ERA5 wind. Typical {std['noun']} morning = medians over {std['noun']} mornings: "
             + ("swell height @ peak period and the direction it comes from." if std.get("ocean") else
                f"total wave height @ {'dominant' if std['peak'] else 'mean'} period and the direction waves come from.")
             + pending,
             fontsize=7.5, color=INK_2, va="top")
    fig.text(0.02, 0.033, f"Coasts are ranked by year-round {std['noun']} %. " + TIDE_NOTE + " " + WIND_NOTE +
             ("" if std.get("ocean") else f" {std['name']} doesn't use wind (lake and sea surf is windswell); "
              f"Ocean good+ does: \u2264 {config.MAX_WIND_KT:g} kt or offshore."),
             fontsize=7.5, color=INK_3, va="top")


def region_page(pdf, region, land):
    """A regional close-up: coasts on a map, the region's surfable season, and a month-by-month table.
    Regions with more than SPLIT_ROWS coasts get the table on a second page; returns the figure(s)."""
    from . import regions
    from .analysis import month_counts

    spec = regions.REGIONS[region]
    tp_ratio = next((r for r in (month_counts(m).tp_ratio for m in range(1, 13)) if r), 1.29)
    std = regions.standard(region)
    rows = sorted(regions.summary(region, tp_ratio), key=lambda r: -r["pct"])  # ranked by year-round surfable %

    if len(rows) <= SPLIT_ROWS:
        fig = _new_page(f"{region}: when and where it's {std['noun']}", _region_subtitle(std, tp_ratio))
        table_h = (len(rows) + 1.3) * 0.022  # the table grows with the coasts; map and chart share the rest
        top = 0.07 + table_h + 0.045
        _region_map(fig, [0.02, top, 0.44, 0.87 - top], spec, rows, land)
        _region_season(fig, [0.53, top + 0.04, 0.44, 0.83 - top - 0.04], rows, std)
        _region_table(fig, [0.02, 0.07, 0.96, table_h], rows, std)
        _region_footer(fig, rows, std)
        pdf.savefig(fig, dpi=220)
        return fig

    overview = _new_page(f"{region}: when and where it's {std['noun']}", _region_subtitle(std, tp_ratio))
    _region_map(overview, [0.02, 0.33, 0.96, 0.55], spec, rows, land)
    _region_season(overview, [0.1, 0.06, 0.8, 0.2], rows, std)
    pdf.savefig(overview, dpi=220)
    table = _new_page(f"{region}: coast by coast", _region_subtitle(std, tp_ratio))
    _region_table(table, [0.02, 0.08, 0.96, 0.8], rows, std)
    _region_footer(table, rows, std)
    pdf.savefig(table, dpi=220)
    return [overview, table]


def build(output, preview_dir=None):
    """Page order: year-round map, the popular-breaks and stretches breakouts, the half-year maps,
    one map per month in config.MONTHS, then the regional pages."""
    from .regions import REGIONS

    land = _land()
    stretches_page = 3

    def maps(title, table, months):
        return lambda pdf: page(pdf, Period(months, title, table), land, stretches_page)

    (year, *halves) = config.SUMMARIES
    pages = [(year[0], maps(*year))]
    pages += [("Popular breaks", spot_page), ("Under-the-radar stretches", stretch_page)]
    pages += [(title, maps(title, table, months)) for title, table, months in halves]
    pages += [(calendar.month_name[m], maps(calendar.month_name[m], f"Best named breaks in {calendar.month_name[m]}", (m,)))
              for m in config.MONTHS]
    pages += [(region, lambda pdf, region=region: region_page(pdf, region, land)) for region in REGIONS]
    assert pages[stretches_page - 1][0] == "Under-the-radar stretches"
    with PdfPages(output, metadata={"Title": "Predictive Surf Atlas"}) as pdf:
        n = 1
        for title, render in pages:
            n += _save(render(pdf), preview_dir, n, title)
    print(f"wrote {output}")


def _save(figs, preview_dir, n, title):
    """Close a rendered page (or several); returns how many pages it was."""
    figs = figs if isinstance(figs, list) else [figs]
    for k, fig in enumerate(figs):
        if preview_dir:
            fig.savefig(f"{preview_dir}/page_{n + k:02d}.png", dpi=80)
        plt.close(fig)
        print(f"  page {n + k}: {title}", flush=True)
    return len(figs)
