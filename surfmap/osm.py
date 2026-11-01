"""Surf spots mapped in OpenStreetMap (sport=surfing), for more places to click in the atlas and
to tell charted coast from uncharted. Data © OpenStreetMap contributors, ODbL 1.0."""

import json
import re

import numpy as np
import requests

from . import config

OVERPASS = "https://overpass-api.de/api/interpreter"
QUERY = '[out:json][timeout:300];(node["sport"="surfing"];way["sport"="surfing"];relation["sport"="surfing"];);out center tags;'
PATH = config.DATA / "osm" / "surfing.json"
NOT_A_SPOT = re.compile(r"shop|store|school|schule|escuela|[ée]cole|scuola|club|rental|hire|hostel|camp|academy|lesson|"
                        r"cent(er|re)|centro|lodge|resort|hotel|house|park|pool", re.I)
SEA_KM = 10  # spots further than this from the sea are river or lake waves


def fetch():
    """One Overpass query for every sport=surfing feature, cached in data/osm/."""
    if PATH.exists():
        return
    r = requests.post(OVERPASS, data={"data": QUERY}, timeout=600,
                      headers={"User-Agent": "surf_tools (https://github.com/jonchaconas/surf_tools)"})
    r.raise_for_status()
    PATH.parent.mkdir(parents=True, exist_ok=True)
    PATH.write_text(r.text)


def spots(land, land_res):
    """[(name, lat, lon)] of named surf spots by the sea: no shops, schools or clubs, and no river
    waves. land is the exposure-ray land raster (0 water), land_res its resolution in degrees."""
    if not PATH.exists():
        return []
    out = []
    reach = int(np.ceil(SEA_KM / 111.2 / land_res))
    h, w = land.shape
    for e in json.loads(PATH.read_text())["elements"]:
        tags = e.get("tags", {})
        name = tags.get("name", "").strip()
        if not name or NOT_A_SPOT.search(name) or any(k in tags for k in ("shop", "amenity", "office", "club", "craft")):
            continue
        lat, lon = (e["lat"], e["lon"]) if "lat" in e else (e["center"]["lat"], e["center"]["lon"])
        r, c = int((90 - lat) / land_res), int((lon + 180) / land_res)
        box = land[max(0, r - reach):min(h, r + reach + 1), [(c + k) % w for k in range(-reach, reach + 1)]]
        if (box == 0).any():
            out.append((name, round(lat, 4), round(lon, 4)))
    return out
