"""The good+ test applied to one morning's sea state."""

import numpy as np

from . import config

G = 9.81


def breaker_height(hs, tp):
    """Komar & Gaughan (1972) breaking height from deep-water Hs (m) and period (s)."""
    return 0.39 * G**0.2 * (tp * hs**2) ** 0.4


def is_good(swell_hs, swell_tp, windsea_hs):
    """Boolean array of good+ mornings. A missing wind-sea partition means no wind sea."""
    windsea_hs = np.nan_to_num(windsea_hs, nan=0.0)
    with np.errstate(invalid="ignore"):
        hb = breaker_height(swell_hs, swell_tp)
        return (
            (((swell_tp > config.MIN_SWELL_PERIOD) & (swell_hs >= config.MIN_SWELL_HEIGHT))
             | ((swell_tp > config.LONG_SWELL_PERIOD) & (swell_hs >= config.LONG_SWELL_HEIGHT)))
            & (hb <= config.MAX_BREAKER)
            & (windsea_hs <= config.MAX_WINDSEA_RATIO * swell_hs)
        )


def has_sea(swell_hs, windsea_hs):
    """Mornings the model has a sea state for. A blank swell partition alongside a wind sea
    means "no swell" (a morning that isn't good), not missing data."""
    return np.isfinite(swell_hs) | np.isfinite(windsea_hs)


def exposed(swell_dir, facing):
    """True where swell arrives from within config.EXPOSURE degrees of the break's facing."""
    diff = np.abs((np.asarray(swell_dir) - facing + 180) % 360 - 180)
    return diff <= config.EXPOSURE


def swell_ok(hs, tp, swell_dir, facing):
    """One swell partition that would make a break good+: big enough for its period, breaking under
    MAX_BREAKER, and arriving within EXPOSURE degrees of the facing."""
    with np.errstate(invalid="ignore"):
        return ((((tp > config.MIN_SWELL_PERIOD) & (hs >= config.MIN_SWELL_HEIGHT))
                 | ((tp > config.LONG_SWELL_PERIOD) & (hs >= config.LONG_SWELL_HEIGHT)))
                & (breaker_height(hs, tp) <= config.MAX_BREAKER) & exposed(swell_dir, facing))


def wind_ok(kt, wind_from, facing, max_kt=config.MAX_WIND_KT):
    """Clean morning wind at a break: light, or blowing from the land (offshore)."""
    with np.errstate(invalid="ignore"):
        offshore = np.abs((np.asarray(wind_from) - facing) % 360 - 180) <= config.OFFSHORE
        return (np.asarray(kt) <= max_kt) | offshore


def break_good(m, facing, max_kt=config.MAX_WIND_KT):
    """Good+ at a break, from points.mornings(): either swell partition passes swell_ok (so a
    groundswell under a local windswell still counts) and the morning wind is clean."""
    return (swell_ok(m["hs1"], m["tp1"], m["dir1"], facing) | swell_ok(m["hs2"], m["tp2"], m["dir2"], facing)) \
        & wind_ok(m["kt"], m["wind_from"], facing, max_kt)


def good_swell(m, facing):
    """(hs, tp, dir) of the partition that made each morning good: the bigger one if both did."""
    ok1 = swell_ok(m["hs1"], m["tp1"], m["dir1"], facing)
    ok2 = swell_ok(m["hs2"], m["tp2"], m["dir2"], facing)
    use2 = ok2 & (~ok1 | (np.nan_to_num(m["hs2"]) > np.nan_to_num(m["hs1"])))
    return tuple(np.where(use2, m[f"{v}2"], m[f"{v}1"]) for v in ("hs", "tp", "dir"))


FT = 3.28084  # feet per metre


def rule_text():
    return (
        f"Good+ morning (Jon's Criteria): swell of at least {config.MIN_SWELL_HEIGHT * FT:.0f} ft ({config.MIN_SWELL_HEIGHT:g} m) "
        f"with peak period over {config.MIN_SWELL_PERIOD:g} s, or {config.LONG_SWELL_HEIGHT * FT:.0f} ft "
        f"({config.LONG_SWELL_HEIGHT:g} m) over {config.LONG_SWELL_PERIOD:g} s; breaking under {config.MAX_BREAKER * FT:.0f} ft "
        f"(Komar–Gaughan from swell height and period), and wind-sea height "
        f"≤ {config.MAX_WINDSEA_RATIO:g}× swell height. "
        f"Judged on one sample near {config.LOCAL_HOUR:02d}:00 local solar time."
    )
