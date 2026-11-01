"""Headless checks for the atlas served by `python -m surfmap gui` on http://127.0.0.1:8765/.

    uv run --with playwright python -I .claude/skills/publish-atlas/check_atlas.py

Screenshots (light, dark, phone) are written next to this file. Prints timings, the worst
mismatch between the click panel and the map over random coastal pixels (must be 0), whether
the phone layout scrolls sideways (must be False), and any console errors (must be clean).
First run: `uv run --with playwright python -m playwright install chromium`.
"""

import asyncio
from pathlib import Path

from playwright.async_api import async_playwright

URL = "http://127.0.0.1:8765/#m=year"  # a link with state: a bare one opens on the forecast near the viewer
HERE = Path(__file__).parent
READY = "document.getElementById('loading').hidden === true"


async def main():
    async with async_playwright() as p:
        browser = await p.chromium.launch()
        errors = []

        def watch(page):
            page.on("pageerror", lambda e: errors.append(f"pageerror: {e}"))
            page.on("console", lambda m: errors.append(f"{m.type}: {m.text}") if m.type in ("error", "warning") else None)

        for scheme in ("light", "dark"):
            page = await browser.new_page(viewport={"width": 1440, "height": 880}, color_scheme=scheme, device_scale_factor=2)
            watch(page)
            await page.goto(URL, wait_until="networkidle")
            await page.wait_for_function(READY, timeout=90000)
            await page.wait_for_timeout(400)
            if scheme == "light":
                print("timing (ms):", await page.evaluate("""() => {
                  let t = performance.now(); recompute(); const all = performance.now() - t;
                  t = performance.now(); renderField(devicePixelRatio); const render = performance.now() - t;
                  return { recompute: Math.round(all), renderField: Math.round(render), points: meta.points.length };
                }"""))
                print("panel vs map:", await page.evaluate("""() => {
                  let worst = 0;
                  for (let n = 0; n < 400; n++) {
                    const c = Math.floor(Math.random() * R.pix.length), q = R.pix[c];
                    const lat = -90 + (Math.floor(q / R.nlon) + 0.5) * R.res, lon = -180 + ((q % R.nlon) + 0.5) * R.res;
                    selected = locate(lon, lat); renderDetail();
                    const shown = parseFloat(document.querySelector('.bigpct').textContent);
                    if (isFinite(F[q])) worst = Math.max(worst, Math.abs(Math.round(F[q]) - shown));
                  }
                  selected = null; renderDetail();
                  return { worst };
                }"""))
            await page.screenshot(path=str(HERE / f"atlas_{scheme}.png"))
            await page.close()

        phone = await browser.new_page(viewport={"width": 400, "height": 860}, is_mobile=True)
        watch(phone)
        await phone.goto(URL, wait_until="networkidle")
        await phone.wait_for_function(READY, timeout=90000)
        print("phone horizontal overflow:", await phone.evaluate("document.documentElement.scrollWidth > window.innerWidth"))
        await phone.screenshot(path=str(HERE / "atlas_phone.png"))
        print("console:", errors or "clean")
        await browser.close()


asyncio.run(main())
