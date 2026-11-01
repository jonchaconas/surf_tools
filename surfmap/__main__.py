"""Download wave model history and render the good+ surf PDF.

    uv run python -m surfmap fetch            # both sources (cached in data/)
    uv run python -m surfmap fetch --source inland
    uv run python -m surfmap build            # writes surf_good_days.pdf
    uv run python -m surfmap discover         # rewrites surfmap/stretches.py
"""

import argparse

from . import config


def main():
    p = argparse.ArgumentParser(prog="surfmap", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download WW3 and Open-Meteo data into data/")
    f.add_argument("--source", choices=("all", "ww3", "inland", "wind", "regions", "copernicus"), default="all")
    b = sub.add_parser("build", help="render one PDF page per month")
    b.add_argument("-o", "--output", default=str(config.ROOT / "surf_good_days.pdf"))
    sub.add_parser("discover", help="rewrite stretches.py: good coast far from the known breaks")
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
    elif args.cmd == "discover":
        from . import discover
        discover.write(discover.find())
    else:
        from . import plot
        plot.build(args.output)


if __name__ == "__main__":
    main()
