# surf_tools

Where is the surf good, and when? This repo turns five years of wave-model history into the
**Predictive Surf Atlas**: a 26-page PDF of world maps for every month, ranked tables of famous
breaks and little-known coast, and close-ups of regions such as the Great Lakes and the
Mediterranean, plus an interactive version where you set your own definition of good.

**[Download the PDF](predictive_surf_atlas.pdf)** · **[Open the interactive atlas](https://jonchaconas.github.io/surf_tools/)**

Everything is computed from open model data. The figures describe the offshore sea state,
not reports from the beach (see [Limitations](#limitations)).

## What's in the PDF

| Pages | Content |
|---|---|
| 1 | **Year-round** world map: share of mornings that are good+, contoured red (rarely) to green (usually), with arrows for the predominant swell direction |
| 2 | **Popular breaks, month by month**: 40 well-known breaks, one per 100 miles of coast |
| 3 | **Under-the-radar stretches, month by month**: coast that is often good but far from the famous breaks |
| 4–5 | **October–March** and **April–September** world maps |
| 6–17 | One world map per month, **November → October** |
| 18–26 | Regional close-ups: Lake Michigan, Lake Superior, the Mediterranean (2 pages), Greece, Caspian Sea, Black Sea, Madagascar, Peru |

Every world map also has:

- the 15 best named breaks for that period, each with its typical morning wind;
- the 7 best under-the-radar stretches, shown as numbered ◆ markers;
- insets for the Great Lakes, the Mediterranean and Black Sea, and the Caspian.

The month-by-month tables show, for each coast:

- the % of good mornings in each month, plus year-round and the worst month;
- the best three months;
- a typical good morning, for example `6 ft @ 13 s SSW, surf ~9 ft`;
- the tidal range;
- the typical morning wind.

Breakout tables are listed by year-round consistency.

## Interactive atlas

[`web/index.html`](web/index.html) is an interactive version of the PDF:

- **Map:** zoom and pan the world map, and switch between year-round, the half years or any
  single month. The contours are drawn from a smooth 0.25° field that runs right up to the
  coastline.
- **Coast that faces the right way:** within about 80 km of land, each 0.25° pixel weights every
  morning by how open it is to that morning's swell direction. Openness is the share of sight
  lines from that direction that run 200 km over water (100 km in lakes and enclosed seas).
  Mainland and large islands block what they face. Islands under 1,000 km², such as the Channel
  Islands, only half-block, since swell wraps around them and through the gaps between them.
- **Click anywhere:** a click shows a month-by-month chart for whatever it lands on:
  - a named break, an under-the-radar stretch or a regional coast, with its typical good
    morning, tide range and morning wind;
  - any stretch of coast, with a compass rose showing how open it is to each direction;
  - any patch of ocean.

  A click inland snaps to the nearest coast.
- **Your own definition of good:** sliders for:
  - minimum swell size and period, and the long-period alternative;
  - the largest surf you want;
  - the morning wind allowed at the breaks;
  - how far off a break's facing the swell may arrive;
  - how much wind chop to allow out at sea.

  The map and the best-spots list update as you go. Presets include small groundswell,
  longboard, shortboard, big days and lake surf.

To make the map work in a browser, each month's mornings for the grid cells are pre-counted into
bins of swell height × period × chop by [`surfmap/export.py`](surfmap/export.py), so any criteria
reduce to sums over bins. For the coastal band, the counts are also split by the compass
direction the swell comes from. The breaks, stretches and regional coasts ship their raw
mornings instead (both swell trains and the wind), so every setting applies to them exactly.

The data in `web/data/` is about 20 MB. It holds:
- a 2° ocean grid and the enclosed-sea points;
- the 0.25° coastal band: for each coastal pixel, which directions are open (found by casting
  rays over a 0.05° land mask) and which model cell it reads;
- every morning at the ~260 marked points.

Run it locally, at 1° resolution, with buttons to download new data and rebuild the PDF:

```sh
uv run python -m surfmap gui
```

`web/` is a static site. [`.github/workflows/pages.yml`](.github/workflows/pages.yml) publishes it to
**https://jonchaconas.github.io/surf_tools/** whenever `web/` changes on `main`. This needs a one-time
setting: Settings → Pages → Source: GitHub Actions.

## What counts as good

Each day is judged on **one model snapshot near 08:00 local solar time**.

**Good+, or Jon's Criteria** (open ocean): the swell is

- at least **4 ft (1.2 m) with a peak period over 9 s**, or at least **3 ft (0.9 m) over 14 s**
  (long-period swell breaks bigger than its height);
- breaking under **16 ft (5 m)** (Komar–Gaughan breaking height);
- the dominant part of the sea: wind-sea height is no more than **0.35×** the swell height.

**At a break** (named breaks, under-the-radar stretches and the regional coasts), the rule is read
from the wave model's ~8 km cell at the break, with real morning wind:

- **either swell train** may qualify, so a small groundswell running under a local windswell
  still counts. Malibu's summer south swells are the example;
- the swell must arrive within **75°** of the direction the break faces;
- the morning wind must be **8 kt or less, or offshore** (from the land, within 45° of straight
  offshore, ERA5). This replaces the wind-sea test, which reads chop far out to sea.

The open-ocean map colours keep the wind-sea test, since they have no shoreline to be offshore of.

**Enclosed seas** (Great Lakes, Mediterranean, Black Sea, Caspian) rarely see groundswell. Their
map contours and inset labels therefore use a relaxed standard: total waves **4–10 ft** at a
**dominant period of at least 6 s**.

The regional pages use the standard that suits each coast:

| Standard | Rule | Used on |
|---|---|---|
| Good lake surf | total waves 4–10 ft at a dominant period ≥ 6 s, aimed at the shore | Lake Michigan, Lake Superior |
| Surfable windswell | total waves 3–10 ft at a mean period ≥ 5 s, aimed at the shore | Mediterranean, Greece, Caspian, Black Sea |
| Good+ | the break rule above | Madagascar, Peru |

All thresholds live in [`surfmap/config.py`](surfmap/config.py) and
[`surfmap/regions.py`](surfmap/regions.py).

## Under-the-radar stretches

[`surfmap/discover.py`](surfmap/discover.py) scans every coastal model cell for coast that meets
all of these:

- good+ with swell aimed at the shore on at least **40%** of mornings year-round, or **65%** across
  its best three months;
- at least **300 km** from every well-known break in [`surfmap/spots.py`](surfmap/spots.py);
- within **150 km** of a town of 5,000 or more.

Where qualifying stretches are less than 200 km apart, only the best is kept. Each stretch is
named after the nearest seaside town and written to
[`surfmap/stretches.py`](surfmap/stretches.py).

## Data sources

| What | Source | Resolution |
|---|---|---|
| Open-ocean waves | NOAA WAVEWATCH III global model, University of Hawaiʻi run, via [PacIOOS ERDDAP](https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global.html) | 1°, plus 0.5° over East and Southeast Asia |
| Breaks, stretches, regional coasts | Météo-France MFWAM via [Copernicus Marine](https://marine.copernicus.eu): the global analysis from Nov 2022 and the WAVERYS reanalysis before, both swell partitions ([`surfmap/points.py`](surfmap/points.py)) | ~8 km (0.2° before Nov 2022) at each break |
| Morning wind at breaks | [ERA5](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-timeseries) 10 m wind, point time series from the Climate Data Store | 0.25°, just offshore of each break |
| Enclosed-sea map grids | MFWAM via the [Open-Meteo Marine API](https://open-meteo.com/en/docs/marine-weather-api) and Copernicus Marine | 0.75–2° grids for the map insets |
| Tides | [EOT20](https://www.seanoe.org/data/00683/79489/) global ocean tide model (DGFI-TUM, CC BY 4.0) | 1/8° |
| Coastlines, towns | [Natural Earth](https://www.naturalearthdata.com) | 1:50m and 1:10m |

**Seasons.** Each month uses its five most recent complete years: November–December use
2021–2025, and January–October use 2022–2026.

**Gaps.** WAVEWATCH III has a few outages, for example 1–15 November 2025. Missing days are
dropped from every source so they all cover the same mornings.

**Peak period.** Where swell peak period isn't available (the Great Lakes, and the
Copernicus-filled years) it is estimated as 1.29× the mean period. That is the median ratio
in Open-Meteo's Mediterranean data.


## Rebuilding it

You need [uv](https://docs.astral.sh/uv/), a free
[Copernicus Marine account](https://data.marine.copernicus.eu/register), and a free
[Climate Data Store account](https://cds.climate.copernicus.eu). For the CDS, accept the ERA5
time-series licence on the dataset's Download tab, then put your personal access token in
`~/.cdsapirc`:

```
url: https://cds.climate.copernicus.eu/api
key: <your token>
```

```sh
uv run copernicusmarine login                       # once; stores your credentials locally

uv run python -m surfmap fetch --source ww3         # open-ocean waves (~2 GB)
uv run python -m surfmap fetch --source points      # waves + ERA5 wind at every break (~1-2 h)
uv run python -m surfmap fetch --source inland      # enclosed-sea grids from Open-Meteo (free-tier quota)
uv run python -m surfmap fetch --source copernicus  # fills anything Open-Meteo hasn't delivered
uv run python -m surfmap discover                   # optional: re-find the under-the-radar stretches
uv run python -m surfmap build                      # writes predictive_surf_atlas.pdf
uv run python -m surfmap export --res 2             # data for the interactive atlas (web/data)
```

**Downloads.** Everything is cached in `data/`, which git ignores, and every fetch can be
re-run: it skips files it already has.

**Open-Meteo limits.** The free tier allows 10,000 weighted calls a day, so the Open-Meteo
fetches throttle themselves. You can skip them entirely: the Copernicus fill covers the same
points from the same model with no quota.

**Tides.** EOT20 isn't fetched automatically. Download the zip from SEANOE and extract
`ocean_tides/{M2,S2,N2,K1,O1}_ocean_eot20.nc` into `data/eot20/ocean_tides/`.

### Code layout

| Module | Role |
|---|---|
| `config.py`, `rating.py` | Seasons, sampling hour, the good+ rule |
| `ww3.py`, `points.py`, `inland.py`, `cmems.py`, `tides.py` | One module per data source |
| `regions.py` | The regional pages' coasts and standards |
| `analysis.py` | Reduces each month to counts, and sums months into periods for the maps and tables |
| `plot.py` | Renders the pages |
| `spots.py`, `stretches.py`, `discover.py` | Named breaks, and the search for little-known coast |

## Claude Code skills

[`.claude/skills/`](.claude/skills) has the project's workflows as skills:

- `rebuild-pdf`: fetch and rebuild the PDF;
- `publish-atlas`: export, browser-test and republish the atlas;
- `add-spot`: add or change a spot;
- `check-spot`: explain a spot's rating month by month.

## Limitations

- **Model cells, not the beach.**
  - Breaks are read at the wave model's ~8 km cell; the map colours come from 1–2° cells.
  - Neither resolves reefs, points, nearshore bathymetry, or the angle at which swell wraps
    into a bay. Breaks that need wrapping swell, such as Lobitos and Cabo Blanco in Peru, are
    likely underrated.
  - The coastal band's island and headland shading is geometric (sight lines), not a wave
    model.
- **Wind at 0.25°.** ERA5 smooths local effects such as sea breezes and canyon winds.
- **One sample a day.** Each morning is judged on a single snapshot near 08:00; a session
  later in the day may differ.
- **Coarse enclosed-sea grids.** The map insets use coarse grids (0.75–2°). The regional
  pages, read at each coast's own ~8 km cell, are more reliable.
