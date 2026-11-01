"""Why a spot rates the way it does: its mornings by month and which part of the rule fails.

    uv run python .claude/skills/check-spot/diagnose.py "Malibu"
    uv run python .claude/skills/check-spot/diagnose.py 34.27 -119.30 215   # lat lon facing
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # the repo root, for surfmap
from surfmap import config, points, rating  # noqa: E402


def find(args):
    if len(args) == 3:
        lat, lon, facing = map(float, args)
        return f"{lat}, {lon}", lat, lon, facing
    name = " ".join(args).lower()
    for _, n, region, lat, lon, facing in points.all_points():
        if n.lower() == name:
            return f"{n} ({region})", lat, lon, facing
    sys.exit(f"no spot named {name!r}; pass lat lon facing instead")


def main():
    label, lat, lon, facing = find(sys.argv[1:])
    m = points.mornings(lat, lon)
    if m is None:
        print("not downloaded yet: fetching waves and wind (~1-2 min)")
        points.OUT.mkdir(parents=True, exist_ok=True)
        points._fetch_waves(lat, lon, facing)
        points._fetch_wind(lat, lon, facing)
        points.mornings.cache_clear()
        m = points.mornings(lat, lon)
    print(f"{label}: facing {facing:g}, wave cell {m['cell'][0]:.3f}, {m['cell'][1]:.3f}")
    print(f"rule: swell >= {config.MIN_SWELL_HEIGHT * rating.FT:.0f} ft over {config.MIN_SWELL_PERIOD:g} s or "
          f">= {config.LONG_SWELL_HEIGHT * rating.FT:.0f} ft over {config.LONG_SWELL_PERIOD:g} s, within "
          f"{config.EXPOSURE} deg of facing; wind <= {config.MAX_WIND_KT:g} kt or offshore\n")

    ft = rating.FT
    print("mon    n   big aimed  wind  GOOD 2ft12s | swell 1              swell 2              wind")
    for month in range(1, 13):
        s = points.season(m, (month,))
        n = len(s["day"])
        if not n:
            print(f"{month:3d}    0")
            continue
        big = aimed = np.zeros(n, bool)
        small = np.zeros(n, bool)
        for k in ("1", "2"):
            hs, tp, d = s["hs" + k], s["tp" + k], s["dir" + k]
            with np.errstate(invalid="ignore"):
                b = (((tp > config.MIN_SWELL_PERIOD) & (hs >= config.MIN_SWELL_HEIGHT))
                     | ((tp > config.LONG_SWELL_PERIOD) & (hs >= config.LONG_SWELL_HEIGHT)))
                small |= (hs * ft >= 2) & (tp >= 12) & rating.exposed(d, facing)
            big = big | b
            aimed = aimed | rating.swell_ok(hs, tp, d, facing)
        wind = rating.wind_ok(s["kt"], s["wind_from"], facing)
        good = rating.break_good(s, facing)
        small &= wind
        pct = lambda a: f"{100 * a.mean():4.0f}%"
        train = lambda k: (f"{np.nanmedian(s['hs' + k]) * ft:4.1f} ft {np.nanmedian(s['tp' + k]):4.1f} s "
                           f"{_compass(s['dir' + k]):>3s}")
        wind_txt = f"{np.median(s['kt']):4.1f} kt {_compass(s['wind_from'])}"
        print(f"{month:3d} {n:4d} {pct(big)} {pct(aimed)} {pct(wind)} {pct(good)} {pct(small)}"
              f" | {train('1')}  {train('2')}  {wind_txt}")


def _compass(deg):
    deg = np.asarray(deg, float)
    deg = deg[np.isfinite(deg)]
    if not len(deg):
        return "-"
    r = np.radians(deg)
    mean = np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean())) % 360
    names = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]
    return names[int((mean + 11.25) // 22.5) % 16]


if __name__ == "__main__":
    main()
