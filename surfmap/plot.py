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
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]

INSETS = [  # title, basins it shows, extent (lon0, lon1, lat0, lat1)
    ("Great Lakes", ("great_lakes",), (-94.2, -75.5, 40.8, 49.3)),
    ("Mediterranean & Black Sea", ("mediterranean", "black_sea"), (-6.5, 42.5, 30.0, 47.5)),
    ("Caspian Sea", ("caspian",), (45.8, 55.5, 36.2, 47.5)),
]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
    "hatch.color": HATCH, "hatch.linewidth": 0.5, "pdf.fonttype": 42,
})


def _land():
    """Natural Earth land with the Great Lakes cut out (the Caspian is already a hole)."""
    land = unary_union(list(shp.Reader(shp.natural_earth("50m", "physical", "land")).geometries()))
    return cfeature.ShapelyFeature([land.difference(inland.geometries()["great_lakes"])], PC)


def _draw(ax, field, nodata, land, rows=slice(None), cols=slice(None)):
    lo, la = LON[cols], LAT[rows]
    cs = ax.contourf(lo, la, field[rows][:, cols], levels=LEVELS, cmap=CMAP, norm=NORM,
                     transform=PC, zorder=1)
    cs.set_rasterized(True)
    nd = ax.contourf(lo, la, nodata[rows][:, cols].astype(float), levels=[0.5, 1.5],
                     colors=[NODATA], hatches=["////"], transform=PC, zorder=1.5)
    nd.set_rasterized(True)
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
    ax.text(0, 0.968, "% of mornings that are good+ with swell inside the break's\n"
            "exposure window; best break per 100 miles of coast",
            fontsize=8, color=INK_2, va="top", linespacing=1.3)
    top = _distinct([s for s in spots if s["basin"] is None], 15)
    y, step = 0.905, 0.0375
    for i, s in enumerate(top, 1):
        ax.text(0.0, y, f"{i}", fontsize=8.5, color=INK_3, va="top")
        ax.text(0.07, y, s["name"], fontsize=8.5, weight="bold", va="top")
        ax.text(0.07, y - 0.017, s["region"], fontsize=7, color=INK_2, va="top")
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
        ax.text(0.09, y - 0.017, s["region"], fontsize=7, color=INK_2, va="top")
        ax.text(1.0, y, f"{s['pct']:.0f}%", fontsize=9, weight="bold", va="top", ha="right")
        y -= step


def _inset(fig, rect, title, basins, extent, field, nodata, land, spots, ready):
    lon0, lon1, lat0, lat1 = extent
    ax = fig.add_axes(rect, projection=ccrs.Mercator(central_longitude=(lon0 + lon1) / 2))
    ax.set_extent(extent, crs=PC)
    pad = 3
    rows = (LAT > lat0 - pad) & (LAT < lat1 + pad)
    cols = (LON > lon0 - pad) & (LON < lon1 + pad)
    _draw(ax, field, nodata, land, rows, cols)
    mine = [s for s in spots if s["basin"] in basins]
    ax.scatter([s["lon"] for s in mine], [s["lat"] for s in mine], s=10, c=INK,
               edgecolors="white", linewidths=0.6, transform=PC, zorder=5)
    for s in mine:
        ax.text(s["lon"], s["lat"] + 0.012 * (lat1 - lat0) * 3, f"{s['name']} {s['pct']:.0f}%",
                fontsize=6.5, ha="center", va="bottom", transform=PC, zorder=6, path_effects=HALO)
    ax.set_title(title, fontsize=9.5, weight="bold", loc="left", pad=4)
    if any(b not in ready for b in basins):
        ax.text(0.5, 0.5, "Open-Meteo data is\nstill downloading", transform=ax.transAxes,
                ha="center", va="center", fontsize=8.5, color=INK_2, zorder=7,
                bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor=COAST, linewidth=0.5))


def _inset_rects(fig, bottom, height, left=0.02, right=0.98):
    """Inset boxes sized to each map's projected aspect, spread with equal gaps across the page."""
    fw, fh = fig.get_size_inches()
    widths = []
    for _, _, (lon0, lon1, lat0, lat1) in INSETS:
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

    ax = fig.add_axes([0.015, 0.375, 0.73, 0.54], projection=ccrs.Robinson(central_longitude=10))
    ax.set_global()
    _draw(ax, field, nodata, land)
    arrows = _swell_arrows(ax, *m.swell_arrows())
    _label_spots(ax, [s for s in spots if s["basin"] is None])
    _stretch_markers(ax, stretches)
    _colorbar(fig, [0.06, 0.345, 0.46, 0.013], ax.collections[0])
    ax.quiverkey(arrows, 0.715, 0.3515, 1, "Predominant swell direction (power-weighted;\n"
                 "omitted where swell comes from mixed directions)", coordinates="figure",
                 labelpos="E", labelsep=0.08, fontproperties=dict(size=7.5), labelcolor=INK_2)

    _table(fig, [0.77, 0.36, 0.215, 0.57], spots, stretches, m.table_title, breakout_page)

    for (title, basins, extent), rect in zip(INSETS, _inset_rects(fig, bottom=0.068, height=0.235)):
        _inset(fig, rect, title, basins, extent, field, nodata, land, spots, m.ready)

    fig.text(0.02, 0.045, rating.rule_text(), fontsize=7.5, color=INK_2, va="top")
    fig.text(0.02, 0.028,
             "Ocean: NOAA WAVEWATCH III global model (Univ. of Hawaiʻi run, PacIOOS ERDDAP), 1° sampling. "
             "Great Lakes, Mediterranean, Black Sea, Caspian: Open-Meteo Marine API (Météo-France MFWAM), "
             "0.75° lakes / 1.5–2° seas. Both judged by the same open-ocean standard.",
             fontsize=7.5, color=INK_2, va="top")
    lakes = (f" Open-Meteo has no swell peak period on the Great Lakes; it is taken as {m.tp_ratio:.2f}× the "
             f"mean swell period, the median ratio in the other seas." if m.tp_ratio else "")
    fig.text(0.02, 0.011,
             "Offshore sea state only: bathymetry, sheltering, tide and wind direction at each break are not modelled. "
             "Break percentages use the nearest seaward model cell within 250 km." + lakes,
             fontsize=7.5, color=INK_3, va="top")
    pdf.savefig(fig, dpi=220)
    return fig


def _compass(deg):
    if np.isnan(deg):
        return "–"
    points = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return f"{points[int((deg + 11.25) // 22.5) % 16]} {deg:.0f}°"


def _ink_on(color):
    """Dark or white text, whichever reads on this fill."""
    r, g, b = (c if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in matplotlib.colors.to_rgb(color))
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.25 else INK


def stretch_page(pdf):
    """One row per discovered stretch: % good+ each month, year-round, and its typical good morning."""
    from . import discover

    rows = Period(tuple(range(1, 13)), "Year-round", "").stretches()
    monthly = {m: {r["num"]: r["pct"] for r in Period((m,), "", "").stretches()} for m in range(1, 13)}

    fig = plt.figure(figsize=(17, 11))
    fig.patch.set_facecolor("white")
    fig.text(0.02, 0.965, "Under-the-radar stretches, month by month", fontsize=24, weight="bold", va="top")
    fig.text(0.02, 0.925,
             f"Coast where mornings are often good+ with swell aimed at the shore (≥{discover.MIN_YEAR}% "
             f"year-round, or ≥{discover.MIN_PEAK}% over its best three months),\nat least {discover.FAR_KM} km "
             f"from the well-known breaks and within {discover.TOWN_KM} km of a town of {discover.TOWN_POP:,}+. "
             f"Numbers match the ◆ markers on the maps.", fontsize=11.5, color=INK_2, va="top",
             linespacing=1.4)

    ax = fig.add_axes([0.02, 0.075, 0.96, 0.8])
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    h = 1 / (len(rows) + 1.5)
    x_month, w_month = 0.335, 0.034
    x_year, x_worst, x_size = x_month + 12 * w_month + 0.012, x_month + 12 * w_month + 0.062, 0.86
    y = 1 - h / 2
    for x, label, ha in [(0.0, "#", "left"), (0.025, "Stretch", "left"), (0.19, "Country", "left"),
                         (x_year + 0.02, "Year", "center"), (x_worst + 0.02, "Worst", "center"),
                         (x_size, "Typical good morning", "left")]:
        ax.text(x, y, label, fontsize=8.5, weight="bold", color=INK_2, ha=ha, va="center")
    for k, mon in enumerate("JFMAMJJASOND"):
        ax.text(x_month + (k + 0.5) * w_month, y, mon, fontsize=8.5, weight="bold", color=INK_2,
                ha="center", va="center")

    for r in rows:
        y -= h
        ax.text(0.0, y, f"◆{r['num']}", fontsize=8.5, color=INK_3, va="center")
        ax.text(0.025, y, r["name"], fontsize=8.5, weight="bold", va="center")
        ax.text(0.19, y, r["region"], fontsize=8, color=INK_2, va="center")
        vals = [monthly[m].get(r["num"], np.nan) for m in range(1, 13)]
        for k, v in enumerate(vals):
            color = CMAP(NORM(v)) if np.isfinite(v) else NODATA
            ax.add_patch(matplotlib.patches.Rectangle((x_month + k * w_month + 0.001, y - h * 0.42),
                                                      w_month - 0.002, h * 0.84, facecolor=color, linewidth=0))
            ax.text(x_month + (k + 0.5) * w_month, y, f"{v:.0f}" if np.isfinite(v) else "–",
                    fontsize=7.5, ha="center", va="center", color=_ink_on(color))
        ax.text(x_year + 0.02, y, f"{r['pct']:.0f}%", fontsize=9, weight="bold", ha="center", va="center")
        ax.text(x_worst + 0.02, y, f"{np.nanmin(vals):.0f}%", fontsize=8.5, ha="center", va="center")
        ax.text(x_size, y, f"{r['hb']:.1f} m, {r['tp']:.0f} s, {_compass(r['dir'])}",
                fontsize=8.5, color=INK_2, va="center")

    fig.text(0.02, 0.05, "Month cells use the map colours (% of mornings good+). Typical good morning = breaker "
             "height (Komar–Gaughan estimate), swell period and the direction swell comes from, as medians "
             "over that stretch's good+ mornings.",
             fontsize=7.5, color=INK_2, va="top")
    fig.text(0.02, 0.033, "Each stretch is read from one WW3 cell 50–120 km offshore. The rule has no wind, "
             "which flatters cold, steady coasts (Chile, Peru, Namibia) and storm coasts (southern Chile, New "
             "Zealand's West Coast); check local wind and access before booking a trip.",
             fontsize=7.5, color=INK_3, va="top")
    pdf.savefig(fig, dpi=220)
    return fig


def build(output, preview_dir=None):
    """Summary pages (config.SUMMARIES), the stretches breakout, then one page per month in config.MONTHS."""
    land = _land()
    breakout = len(config.SUMMARIES) + 1
    periods = [(title, table, months) for title, table, months in config.SUMMARIES]
    periods += [(calendar.month_name[m], f"Best named breaks in {calendar.month_name[m]}", (m,))
                for m in config.MONTHS]
    with PdfPages(output) as pdf:
        for i, (title, table, months) in enumerate(periods, 1):
            if i == breakout:
                _save(stretch_page(pdf), preview_dir, i, "Under-the-radar stretches")
            n = i + (i >= breakout)
            _save(page(pdf, Period(months, title, table), land, breakout), preview_dir, n, title)
    print(f"wrote {output}")


def _save(fig, preview_dir, n, title):
    if preview_dir:
        fig.savefig(f"{preview_dir}/page_{n:02d}.png", dpi=80)
    plt.close(fig)
    print(f"  page {n}: {title}", flush=True)
