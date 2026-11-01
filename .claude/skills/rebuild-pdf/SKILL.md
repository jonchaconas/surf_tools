---
name: rebuild-pdf
description: Download any missing wave/wind data and rebuild predictive_surf_atlas.pdf, then check the pages render. Use when asked to update, regenerate, refresh or send the surf PDF, or after changing the rule, spots or regions.
---

# Rebuild the surf PDF

The PDF is `predictive_surf_atlas.pdf` at the repo root (26 pages: year-round, popular breaks,
under-the-radar stretches, two half-years, twelve months, regional pages).

## 1. Data

Everything is cached under `data/` and every fetch skips what it already has, so re-running is
cheap. Only fetch what the change needs:

```sh
uv run python -m surfmap fetch --source ww3         # open-ocean waves (PacIOOS; slow server, be patient)
uv run python -m surfmap fetch --source points      # breaks, stretches, regional coasts: Copernicus + ERA5
uv run python -m surfmap fetch --source copernicus  # enclosed-sea grids Open-Meteo hasn't delivered
```

- **Copernicus Marine** needs `uv run copernicusmarine login` (done once by the user).
- **ERA5** needs `~/.cdsapirc`. The user writes the key; never ask them to paste it in chat,
  and never print the file.
- **Open-Meteo** (`--source regions`, `--source inland`) has a free-tier quota. Never route
  around it (VPNs, extra IPs). Copernicus covers the same points without a quota.
- A stalled WW3 fetch is usually PacIOOS being slow, not a bug. Check the log before changing
  code.

`points` is the long one (about an hour for all ~260 points). Run it in the background with a
log, for example `nohup uv run python -c "from surfmap import points; points.fetch(workers=6)" >
data/fetch_points.log 2>&1 &`. The HDF5-DIAG noise in that log is harmless; real failures print
as `... failed:` and the final line is `done: N, failed: M`.

## 2. Build

```sh
uv run python -m surfmap build        # writes predictive_surf_atlas.pdf (several minutes)
```

## 3. Check before handing it over

Render the pages to PNG and look at them. Don't trust the build exiting 0:

```sh
uv run --with pymupdf python - <<'EOF'
import fitz  # PyMuPDF
doc = fitz.open("predictive_surf_atlas.pdf")
for i, page in enumerate(doc, 1):
    page.get_pixmap(dpi=60).save(f"/tmp/page_{i:02d}.png")
print(len(doc), "pages")
EOF
```

Look for:

- maps with no red anywhere, which means a contour reprojection failed;
- overlapping labels or text running off the page;
- regional pages showing "still downloading" dashes;
- table rows missing for breaks whose `points` download failed.

Then commit the PDF with the code change and push (the repo is public:
github.com/jonchaconas/surf_tools).

## Conventions the user has set

- Red→green contours.
- Gaza and Beirut, never Tel Aviv. Montara, not Mavericks.
- Breakout tables are sorted by year-round consistency.
- Swell in ft, wind in knots, tide ranges shown.
