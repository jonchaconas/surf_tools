"""Download wave model history and render the good+ surf PDF.

    uv run python -m surfmap fetch            # both sources (cached in data/)
    uv run python -m surfmap fetch --source inland
    uv run python -m surfmap fetch --source points   # Copernicus waves + ERA5 wind at each break
    uv run python -m surfmap fetch --source cdip     # California nearshore waves (CDIP MOP)
    uv run python -m surfmap build            # writes predictive_surf_atlas.pdf
    uv run python -m surfmap discover         # rewrites surfmap/stretches.py
    uv run python -m surfmap gui              # the interactive atlas, locally
"""

import argparse

from . import config


def main():
    p = argparse.ArgumentParser(prog="surfmap", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download WW3 and Open-Meteo data into data/")
    f.add_argument("--source", choices=("all", "ww3", "inland", "wind", "regions", "copernicus", "points", "cdip", "osm"), default="all")
    b = sub.add_parser("build", help="render one PDF page per month")
    b.add_argument("-o", "--output", default=str(config.ROOT / "predictive_surf_atlas.pdf"))
    sub.add_parser("discover", help="rewrite stretches.py: good coast far from the known breaks")
    e = sub.add_parser("export", help="write web/data/ for the interactive atlas")
    e.add_argument("--res", type=int, default=2, help="ocean grid in degrees: 2 for the shareable page, 1 locally")
    g = sub.add_parser("gui", help="run the interactive atlas locally")
    g.add_argument("--port", type=int, default=8765)
    g.add_argument("--res", type=int, default=1, help="grid used if the app data has to be exported first")
    g.add_argument("--no-browser", action="store_true")
    args = p.parse_args()

    if args.cmd == "fetch":
        from . import inland, regions, wind, ww3
        if args.source in ("all", "ww3"):
            ww3.fetch()
        if args.source in ("all", "wind"):
            wind.fetch()
        if args.source in ("all", "regions"):
            regions.fetch()
        if args.source in ("all", "inland"):
            inland.fetch()
        if args.source == "copernicus":  # fills whatever Open-Meteo hasn't delivered; needs a login
            from . import cmems
            cmems.fill()
        if args.source in ("all", "points"):  # needs a Copernicus login and a CDS key
            from . import points
            points.fetch()
        if args.source in ("all", "cdip"):  # California nearshore waves, no login
            from . import cdip
            cdip.fetch()
        if args.source in ("all", "osm"):  # OpenStreetMap surf spots, for the atlas
            from . import osm
            osm.fetch()
    elif args.cmd == "discover":
        from . import discover
        discover.write(discover.find())
    elif args.cmd == "export":
        from . import export
        export.export(args.res)
    elif args.cmd == "gui":
        from . import gui
        gui.serve(args.port, args.res, not args.no_browser)
    else:
        from . import plot
        plot.build(args.output)


if __name__ == "__main__":
    main()
