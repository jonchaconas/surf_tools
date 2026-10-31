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
from .analysis import LAT, LON, Month

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
    ax.add_feature(land, facecolor=LAND, edgecolor="none", zorder=2)
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), edgecolor=COAST, linewidth=0.3, zorder=3)
    ax.add_feature(cfeature.LAKES.with_scale("50m"), facecolor="none", edgecolor=COAST,
                   linewidth=0.3, zorder=3)
    ax.spines["geo"].set_edgecolor(COAST)
    ax.spines["geo"].set_linewidth(0.5)


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


def _table(fig, rect, spots, month_name):
    ax = fig.add_axes(rect)
    ax.axis("off")
    ax.text(0, 1, f"Best named breaks in {month_name}", fontsize=12, weight="bold", va="top")
    ax.text(0, 0.968, "% of mornings that are good+ with swell\ninside the break's exposure window",
            fontsize=8, color=INK_2, va="top", linespacing=1.3)
    top = [s for s in spots if s["basin"] is None][:24]
    y, step = 0.905, 0.0375
    for i, s in enumerate(top, 1):
        ax.text(0.0, y, f"{i}", fontsize=8.5, color=INK_3, va="top")
        ax.text(0.07, y, s["name"], fontsize=8.5, weight="bold", va="top")
        ax.text(0.07, y - 0.017, s["region"], fontsize=7, color=INK_2, va="top")
        ax.text(1.0, y, f"{s['pct']:.0f}%", fontsize=9, weight="bold", va="top", ha="right")
        y -= step


def _inset(fig, rect, title, basins, extent, field, nodata, land, spots):
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
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", spacing="uniform")
    cb.outline.set_visible(False)
    cb.set_ticks(LEVELS, labels=[f"{v}%" for v in LEVELS])
    cb.ax.tick_params(labelsize=8, length=0, colors=INK_2)
    cax.set_title("Share of mornings with good-or-better surf", fontsize=8.5, color=INK_2, loc="left", pad=3)
    fig.patches.append(matplotlib.patches.Rectangle(
        (rect[0] + rect[2] + 0.02, rect[1]), 0.022, rect[3], transform=fig.transFigure,
        facecolor=NODATA, hatch="////", edgecolor=HATCH, linewidth=0))
    fig.text(rect[0] + rect[2] + 0.046, rect[1] + rect[3] / 2, "No model data (sea ice, or a\nsea neither source covers)",
             fontsize=7.5, color=INK_2, va="center", linespacing=1.2)


def _years_text(m):
    days = m.days
    years = sorted({d.year for d in days})
    full = calendar.monthrange(years[-1], m.month)[1] * len(years)
    span = f"{years[0]}–{years[-1]}"
    gap = "" if len(days) == full else f"; {full - len(days)} days missing from the archive"
    return f"{len(years)} seasons, {span} ({len(days)} mornings{gap})"


def page(pdf, month, land):
    m = Month(month)
    field, nodata = m.field()
    spots = m.spots()
    name = calendar.month_name[month]

    fig = plt.figure(figsize=(17, 11))
    fig.patch.set_facecolor("white")
    fig.text(0.02, 0.965, f"{name}: where the surf is good", fontsize=24, weight="bold", va="top")
    fig.text(0.02, 0.925, f"Share of mornings with good-or-better surf in the wave-model record, "
             f"{_years_text(m)}", fontsize=11.5, color=INK_2, va="top")

    ax = fig.add_axes([0.015, 0.375, 0.73, 0.54], projection=ccrs.Robinson(central_longitude=10))
    ax.set_global()
    _draw(ax, field, nodata, land)
    _label_spots(ax, [s for s in spots if s["basin"] is None])
    _colorbar(fig, [0.06, 0.345, 0.46, 0.013], ax.collections[0])

    _table(fig, [0.77, 0.36, 0.215, 0.57], spots, name)

    for (title, basins, extent), rect in zip(INSETS, _inset_rects(fig, bottom=0.068, height=0.235)):
        _inset(fig, rect, title, basins, extent, field, nodata, land, spots)

    fig.text(0.02, 0.045, rating.rule_text(), fontsize=7.5, color=INK_2, va="top")
    fig.text(0.02, 0.028,
             "Ocean: NOAA WAVEWATCH III global model (Univ. of Hawaiʻi run, PacIOOS ERDDAP), 1° sampling. "
             "Great Lakes, Mediterranean, Black Sea, Caspian: Open-Meteo Marine API (Météo-France MFWAM), "
             "0.75° lakes / 1.5–2° seas. Both judged by the same open-ocean standard.",
             fontsize=7.5, color=INK_2, va="top")
    ratio = next((b["tp_ratio"] for b in m.basins.values() if b["tp_ratio"]), None)
    lakes = (f" Open-Meteo has no swell peak period on the Great Lakes; it is taken as {ratio:.2f}× the mean "
             f"swell period, the median ratio in the other seas." if ratio else "")
    fig.text(0.02, 0.011,
             "Offshore sea state only: bathymetry, sheltering, tide and wind direction at each break are not modelled. "
             "Break percentages use the nearest seaward model cell within 250 km." + lakes,
             fontsize=7.5, color=INK_3, va="top")
    pdf.savefig(fig, dpi=220)
    return fig


def build(output, preview_dir=None):
    land = _land()
    with PdfPages(output) as pdf:
        for month in config.MONTHS:
            fig = page(pdf, month, land)
            if preview_dir:
                fig.savefig(f"{preview_dir}/page_{month:02d}.png", dpi=80)
            plt.close(fig)
            print(f"  page {calendar.month_name[month]}", flush=True)
    print(f"wrote {output}")
