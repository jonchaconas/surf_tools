# surf-forecast-tools

Maps of how often the surf is good around the world, built from wave-model history.
`surf_good_days.pdf` has one page per month. Each page has:

- a world map with filled contours of the share of mornings that were good-or-better (good+)
- insets for the Great Lakes, the Mediterranean and Black Sea, and the Caspian
- a ranked table of named breaks

## What counts as good+

Each day is judged on one model snapshot near 08:00 local solar time. A morning is
good+ when the swell:

- has a peak period of at least 10 s,
- breaks between 1.5 and 5 m (Komar–Gaughan breaking height from swell height and period),
- dominates the sea state, with wind-sea height no more than 0.35× the swell height.

The thresholds live in [`surfmap/config.py`](surfmap/config.py). For the break table,
the swell must also arrive within 75° of the direction the break faces.

## Data

| Area | Source | Sampling |
|---|---|---|
| Open ocean | NOAA WAVEWATCH III global (University of Hawaiʻi run, [PacIOOS ERDDAP](https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global.html)) | 1° grid |
| Great Lakes, Mediterranean, Black Sea, Caspian | [Open-Meteo Marine API](https://open-meteo.com/en/docs/marine-weather-api) (Météo-France MFWAM) | 0.75° on the lakes, 1.5–2° on the seas |

Open-Meteo's swell fields start in October 2021, so each month uses its five most
recent complete seasons. WW3 is missing 1–15 November 2025; those days are dropped
from both sources. Open-Meteo has no swell peak period on the Great Lakes, so there
it is estimated from the mean swell period, using the median peak/mean ratio from the
other seas.

## Usage

Requires [uv](https://docs.astral.sh/uv/).

```sh
uv run python -m surfmap fetch    # downloads into data/ (cached, resumable)
uv run python -m surfmap build    # writes surf_good_days.pdf
```

The WW3 download is about 280 MB for two months. Open-Meteo's free tier limits weighted
calls to 5,000 per hour and 10,000 per day. The inland fetch throttles itself and can be
re-run safely: it skips files it has already downloaded.

## Limitations

This measures the offshore sea state only. Local bathymetry, swell shadowing, tide, and
wind direction at the break are not modelled. Treat the break table as a measure of
swell exposure, not a surf report.
