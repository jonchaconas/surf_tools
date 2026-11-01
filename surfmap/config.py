"""Seasons, sampling time and the good+ rule, shared by every data source."""

from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

MONTHS = (11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)  # page order: the year ahead from November

# Summary pages: year-round first, then the half years. (page title, table title, months).
SUMMARIES = (
    ("Year-round", "Best named breaks year-round", tuple(range(1, 13))),
    ("October–March", "Best named breaks, October–March", (10, 11, 12, 1, 2, 3)),
    ("April–September", "Best named breaks, April–September", (4, 5, 6, 7, 8, 9)),
)
N_SEASONS = 5  # Open-Meteo only has swell partitions from Oct 2021, which caps the record


def seasons(month):
    """The N_SEASONS most recent years in which this calendar month is complete."""
    today = date.today()
    last = today.year if month < today.month else today.year - 1
    return range(last - N_SEASONS + 1, last + 1)

# Each day is judged on one snapshot near this local solar hour (dawn-patrol conditions).
LOCAL_HOUR = 8

# A morning is good+ when all of these hold.
MIN_SWELL_PERIOD = 10.0   # s, swell peak period
MIN_BREAKER = 1.5         # m, Komar-Gaughan breaking height of the swell
MAX_BREAKER = 5.0         # m, above this most breaks close out or need a big-wave crew
MAX_WINDSEA_RATIO = 0.35  # wind-sea Hs / swell Hs; low means a clean, swell-dominated sea

# Spot rankings also require the swell to come from within this many degrees of the facing.
EXPOSURE = 75
