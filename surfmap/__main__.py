"""Download wave model history and render the good+ surf PDF.

    uv run python -m surfmap fetch            # both sources (cached in data/)
    uv run python -m surfmap fetch --source inland
    uv run python -m surfmap build            # writes surf_good_days.pdf
"""

import argparse

from . import config


def main():
    p = argparse.ArgumentParser(prog="surfmap", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download WW3 and Open-Meteo data into data/")
    f.add_argument("--source", choices=("all", "ww3", "inland"), default="all")
    b = sub.add_parser("build", help="render one PDF page per month")
    b.add_argument("-o", "--output", default=str(config.ROOT / "surf_good_days.pdf"))
    args = p.parse_args()

    if args.cmd == "fetch":
        from . import inland, ww3
        if args.source in ("all", "ww3"):
            ww3.fetch()
        if args.source in ("all", "inland"):
            inland.fetch()
    else:
        from . import plot
        plot.build(args.output)


if __name__ == "__main__":
    main()
