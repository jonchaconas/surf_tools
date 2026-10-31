"""Seasons, sampling time and the good+ rule, shared by every data source."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# Open-Meteo only has swell partitions from Oct 2021, so every source uses these seasons.
YEARS = range(2021, 2026)
MONTHS = (11, 12)

# Each day is judged on one snapshot near this local solar hour (dawn-patrol conditions).
LOCAL_HOUR = 8

# A morning is good+ when all of these hold.
MIN_SWELL_PERIOD = 10.0   # s, swell peak period
MIN_BREAKER = 1.5         # m, Komar-Gaughan breaking height of the swell
MAX_BREAKER = 5.0         # m, above this most breaks close out or need a big-wave crew
MAX_WINDSEA_RATIO = 0.35  # wind-sea Hs / swell Hs; low means a clean, swell-dominated sea

# Spot rankings also require the swell to come from within this many degrees of the facing.
EXPOSURE = 75
