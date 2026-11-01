"""Live wave buoys for the atlas map: NOAA NDBC's latest observations, refreshed hourly by
.github/workflows/buoys.yml (buoy data can't be fetched from the page itself: NDBC sends no CORS headers).

Writes web/data/buoys_live.json:
    {"generated": ISO time, "buoys": [{"id", "name", "lat", "lon", "t": ISO UTC of the reading,
     "hs", "dpd", "apd", "mwd": all waves (m, dominant and average period s, direction °),
     "swh", "swp", "swd", "wwh", "wwp", "wwd": swell and wind waves from the spectral file (m, s, compass) or null,
     "kt", "wdir", "gust": wind (kt, °), "atmp", "wtmp": air and water (°C),
     "hist": [[hours ago, hs, swh, swp], ...] over the last day}]}
Only buoys reporting waves within MAX_AGE_H hours are kept. Data: NOAA National Data Buoy Center.
"""

import json
import re
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import requests

from . import config

LATEST = "https://www.ndbc.noaa.gov/data/latest_obs/latest_obs.txt"
SPEC = "https://www.ndbc.noaa.gov/data/realtime2/{}.spec"
STATIONS = "https://www.ndbc.noaa.gov/activestations.xml"
OUT = config.ROOT / "web" / "data" / "buoys_live.json"
MAX_AGE_H = 6
SPEC_BYTES = 3000  # the newest ~36 hourly rows of a station's spectral summary
KT = 1.943844


def _num(x):
    try:
        v = float(x)
    except ValueError:  # MM, N/A
        return None
    return None if v in (-99, 99, 999, 9999) else v


def _time(fields):
    y, mo, d, h, mi = (int(float(x)) for x in fields)
    return datetime(y, mo, d, h, mi, tzinfo=timezone.utc)


def _spectral(s, station, now):
    """(latest swell and wind-wave split, last day's [hours ago, hs, swh, swp]) from the spectral file's head."""
    r = s.get(SPEC.format(station), headers={"Range": f"bytes=0-{SPEC_BYTES}"}, timeout=30)
    if r.status_code not in (200, 206):
        return {}, []
    rows = [line.split() for line in r.text.splitlines() if not line.startswith("#")]
    rows = [f for f in rows if len(f) == 15]  # the cut-off last line drops out
    latest, hist = {}, []
    for f in rows:
        ago = (now - _time(f[:5])).total_seconds() / 3600
        if ago > 24:
            break
        if not latest:
            latest = dict(swh=_num(f[6]), swp=_num(f[7]), wwh=_num(f[8]), wwp=_num(f[9]),
                          swd=None if f[10] in ("MM", "N/A") else f[10], wwd=None if f[11] in ("MM", "N/A") else f[11])
        hist.append([round(ago, 1), _num(f[5]), _num(f[6]), _num(f[7])])
    return latest, hist


def live():
    s = requests.Session()
    s.headers["User-Agent"] = "surf_tools (https://github.com/jonchaconas/surf_tools)"
    names = dict(re.findall(r'<station id="([^"]+)"[^>]*?name="([^"]*)"', s.get(STATIONS, timeout=60).text))
    lines = s.get(LATEST, timeout=60).text.splitlines()
    head = lines[0].lstrip("#").split()
    now = datetime.now(timezone.utc)
    buoys = []
    for line in lines[2:]:
        f = dict(zip(head, line.split()))
        if len(f) < len(head) or _num(f["WVHT"]) is None:
            continue
        t = _time([f["YYYY"], f["MM"], f["DD"], f["hh"], f["mm"]])
        if (now - t).total_seconds() > MAX_AGE_H * 3600:
            continue
        kt = lambda v: None if v is None else round(v * KT, 1)
        buoys.append(dict(id=f["STN"], name=names.get(f["STN"], ""), lat=float(f["LAT"]), lon=float(f["LON"]),
                          t=t.strftime("%Y-%m-%dT%H:%MZ"), hs=_num(f["WVHT"]), dpd=_num(f["DPD"]), apd=_num(f["APD"]),
                          mwd=_num(f["MWD"]), kt=kt(_num(f["WSPD"])), wdir=_num(f["WDIR"]), gust=kt(_num(f["GST"])),
                          atmp=_num(f["ATMP"]), wtmp=_num(f["WTMP"])))
    with ThreadPoolExecutor(8) as pool:
        for b, (latest, hist) in zip(buoys, pool.map(lambda b: _safe_spectral(s, b["id"], now), buoys)):
            b.update(dict(swh=None, swp=None, swd=None, wwh=None, wwp=None, wwd=None) | latest, hist=hist)
    OUT.write_text(json.dumps(dict(generated=now.strftime("%Y-%m-%dT%H:%MZ"), buoys=buoys), separators=(",", ":")))
    print(f"live buoys: {len(buoys)} reporting waves in the last {MAX_AGE_H} h; wrote {OUT}", flush=True)


def _safe_spectral(s, station, now):
    try:
        return _spectral(s, station, now)
    except (requests.RequestException, ValueError):
        return {}, []
