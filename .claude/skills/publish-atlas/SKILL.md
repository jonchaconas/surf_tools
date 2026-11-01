---
name: publish-atlas
description: Export web/data for the interactive Predictive Surf Atlas, test the page in a headless browser, republish it to its claude.ai artifact URL and push (which also deploys GitHub Pages). Use after changing web/index.html, surfmap/export.py, the rule, or any data the atlas shows.
---

# Export, test and publish the atlas

The atlas is `web/index.html` plus `web/data/`. It is published in two places:

- a private claude.ai artifact, https://claude.ai/artifact/5Wa8fpAjfn87AZdXgummjP — keep this
  URL; never publish a new artifact for it;
- GitHub Pages, https://jonchaconas.github.io/surf_tools/, deployed by
  `.github/workflows/pages.yml` on every push that touches `web/`.

## 1. Export

```sh
uv run python -m surfmap export --res 2    # ~35 s; prints each file's size
```

Lines like `not downloaded yet, left out: <name>` mean that point's `fetch --source points`
download is missing. Fetch it first, or say so when you publish.

Each file must stay under 16 MB; `grid.b64.txt` is the largest at ~12 MB.

## 2. Test in a browser

Start the local server in the background, then run the checks:

```sh
uv run python -m surfmap gui --no-browser          # serves http://127.0.0.1:8765/
uv run --with playwright python -I .claude/skills/publish-atlas/check_atlas.py
```

`check_atlas.py` loads the page in light, dark and phone layouts and writes screenshots next to
itself. It also checks that the click panel matches the map at 400 random coastal pixels, and
prints recompute and render timings. It must print `console: clean` and `worst: 0`.

**Look at the screenshots.** Contours should run smoothly to the coast, the map should be
cropped to 79°N–68.5°S with the legend underneath, and nothing should overflow at phone width.
Stop the server afterwards.

## 3. Publish

Publish `web/index.html` with the Artifact tool, passing the URL above and the data files as
`files`:

```
data/meta.json, data/land.json, data/grid.b64.txt, data/points.b64.txt,
data/raster.json, data/coast.b64.txt
```

Every path maps to `web/<same path>`. Files you leave out keep their last published version,
so a page-only change needs only `file_path`.

## 4. Commit and push

Commit `web/` together with the code that produced it, then push. If the Pages workflow fails
with "Ensure GitHub Pages has been enabled", the user must set Settings → Pages → Source:
GitHub Actions. The workflow itself is fine.

## Page constraints (artifact host)

- No `<html>`, `<head>` or `<body>` in `web/index.html`. The host adds the skeleton;
  `surfmap/gui.py` and the Pages workflow add the same one locally.
- Scripts only from the CDN allowlist (d3 7.8.5 from cdnjs). Binary data ships as base64
  gzip in `.txt`, because `.gz` isn't served.
