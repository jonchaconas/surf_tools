---
name: check-spot
description: Explain why a surf spot rates the way it does (e.g. "why is Malibu never good", "Ventura should be better in winter"): break its mornings down by month and by which part of the good+ rule fails, and compare criteria. Use when a spot's numbers look wrong or the user questions a rating.
---

# Diagnose a spot's rating

Run the breakdown for any point with data in `data/points/` (named breaks, stretches, regional
coasts), or for any coordinates once fetched:

```sh
uv run python .claude/skills/check-spot/diagnose.py "Malibu"
uv run python .claude/skills/check-spot/diagnose.py 34.27 -119.30 215   # lat lon facing (fetches if needed)
```

It prints, per month:

| Column | Meaning |
|---|---|
| `n` | mornings in the season, with both waves and wind |
| `big` | either swell train is big enough for its period (any direction) |
| `aimed` | … and arrives within the exposure window |
| `wind` | morning wind is light (≤ `MAX_WIND_KT`) or offshore |
| `GOOD` | the good+ rule at the break |
| `2ft12s` | the "small groundswell" preset (2 ft @ 12 s), for spots that live on small long-period swell |

It also prints the typical swell trains (height, period and direction of swell 1 and swell 2)
and the morning wind, so you can see which train the spot depends on.

## Reading it

- **`big` high but `aimed` low:** the swell comes from outside the window. Check `facing`, or
  the spot is sheltered from the dominant swell; Malibu blocks the winter W/NW.
- **`aimed` fine but `GOOD` low:** the wind is the limit. Compare with the user's local
  knowledge; offshore wind always counts as clean.
- **Low everywhere, but locals say it's good:** the spot probably works on small long-period
  swell. Check the `2ft12s` column, and suggest the "Small groundswell" preset rather than
  changing the global rule.
- The model cell is ~8 km, so reefs, points, bathymetry and refraction into bays aren't
  resolved. Say so instead of tuning numbers to match one spot.

User feedback so far (Nov 2026), to check against after rule changes:

- Ventura is good in winter: good+ about 30% of mornings Dec–Feb.
- Malibu works on summer south swells: about 30% Jun–Jul with the "Small groundswell" preset.
