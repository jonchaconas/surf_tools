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
            (swell_tp >= config.MIN_SWELL_PERIOD)
            & (hb >= config.MIN_BREAKER)
            & (hb <= config.MAX_BREAKER)
            & (windsea_hs <= config.MAX_WINDSEA_RATIO * swell_hs)
        )


def exposed(swell_dir, facing):
    """True where swell arrives from within config.EXPOSURE degrees of the break's facing."""
    diff = np.abs((np.asarray(swell_dir) - facing + 180) % 360 - 180)
    return diff <= config.EXPOSURE


def rule_text():
    return (
        f"Good+ morning: swell peak period ≥ {config.MIN_SWELL_PERIOD:g} s, "
        f"breaking height {config.MIN_BREAKER:g}–{config.MAX_BREAKER:g} m "
        f"(Komar–Gaughan, from swell height and period), and wind-sea height "
        f"≤ {config.MAX_WINDSEA_RATIO:g}× swell height. "
        f"Judged on one sample near {config.LOCAL_HOUR:02d}:00 local solar time."
    )
