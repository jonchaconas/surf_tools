# surf-forecast-tools

Maps of how often the surf is good around the world, built from wave-model history.
`surf_good_days.pdf` opens with three summary pages (October–March, April–September,
year-round) and a month-by-month breakout of under-the-radar stretches. It then has one
page per month, November through October. Each map page has:

- a world map with filled contours of the share of mornings that were good-or-better
  (good+), from red (rarely) to green (usually), plus arrows showing the predominant
  swell direction
- insets for the Great Lakes, the Mediterranean and Black Sea, and the Caspian
- a ranked table of named breaks (the best one per 100 miles of coast)
- numbered ◆ markers and a short table for under-the-radar stretches

## Under-the-radar stretches

[`surfmap/discover.py`](surfmap/discover.py) scans every coastal model cell for stretches
that meet all of these:

- often good+ with swell aimed at the shore: at least 40% of mornings year-round, or
  65% across their best three months;
- at least 300 km from every well-known break in [`surfmap/spots.py`](surfmap/spots.py);
- within 150 km of a town of 5,000 or more.

Where qualifying stretches are less than 200 km apart, only the best is kept. The result
is written to [`surfmap/stretches.py`](surfmap/stretches.py).

## What counts as good+

Each day is judged on one model snapshot near 08:00 local solar time. A morning is
good+ when the swell:

- has a peak period of at least 10 s,
- breaks between 1.5 and 5 m (Komar–Gaughan breaking height from swell height and period),
- dominates the sea state, with wind-sea height no more than 0.35× the swell height.

The thresholds live in [`surfmap/config.py`](surfmap/config.py). For the break table,
the swell must also arrive within 75° of the direction the break faces.

Each arrow is the average of a cell's morning swell directions, weighted by wave power
(height² × period). Cells with no clear predominant direction get no arrow: that is,
where the weighted mean vector is under 30% of the total weight.

## Data

| Area | Source | Sampling |
|---|---|---|
| Open ocean | NOAA WAVEWATCH III global (University of Hawaiʻi run, [PacIOOS ERDDAP](https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global.html)) | 1° grid |
| Great Lakes, Mediterranean, Black Sea, Caspian | [Open-Meteo Marine API](https://open-meteo.com/en/docs/marine-weather-api) (Météo-France MFWAM) | 0.75° on the lakes, 1.5–2° on the seas |

Open-Meteo's swell fields start in October 2021, so each month uses its five most
recent complete seasons: 2021–2025 for October to December, 2022–2026 for January to
September. WW3 has a few outages, including 1–15 November 2025, 1 day of July 2022 and
7 days of October 2025. Days missing from WW3 are dropped from both sources. Open-Meteo has no swell peak period on the Great Lakes, so there
it is estimated from the mean swell period, using the median peak/mean ratio from the
other seas.

## Usage

Requires [uv](https://docs.astral.sh/uv/).

```sh
uv run python -m surfmap fetch    # downloads into data/ (cached, resumable)
uv run python -m surfmap build    # writes surf_good_days.pdf
uv run python -m surfmap discover # rewrites surfmap/stretches.py from the downloaded data
```

The WW3 download is about 1.7 GB for all twelve months.

The enclosed seas need about 20,000 Open-Meteo weighted calls. The free tier allows 5,000
per hour and 10,000 per day, so this fetch takes about two days. It fills one month at a
time, waits out throttling, and can be re-run safely: it skips files it has already
downloaded. A month whose enclosed-sea data is incomplete still renders, but its insets
say the data is still downloading.

## Limitations

This measures the offshore sea state only. Local bathymetry, swell shadowing, tide, and
wind direction at the break are not modelled. Treat the break table as a measure of
swell exposure, not a surf report.
