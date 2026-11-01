"""Under-the-radar surf coast found by discover.py (regenerate with `python -m surfmap discover`).

Each entry is (name, region, lat, lon, facing, cell): a coast point, the compass direction
it faces out to sea, and the WW3 1 degree cell its figures come from. Numbering on the
maps follows this order (best year-round first).
"""

STRETCHES = [
    ('Ancud coast', 'Chile', -41.45, -73.86, 244, (-41.5, -74.0)),  # 77% year-round, 87% peak season
    ('Antofagasta coast', 'Chile', -24.33, -70.55, 248, (-24.5, -71.0)),  # 77% year-round, 83% peak season
    ('La Ligua coast', 'Chile', -31.5, -71.58, 269, (-31.5, -72.0)),  # 76% year-round, 80% peak season
    ('Camaná coast', 'Peru', -16.61, -72.79, 193, (-17.5, -73.0)),  # 74% year-round, 80% peak season
    ('Toltén coast', 'Chile', -38.51, -73.52, 271, (-38.5, -74.0)),  # 72% year-round, 76% peak season
    ('La Serena coast', 'Chile', -29.35, -71.35, 255, (-29.5, -72.0)),  # 64% year-round, 71% peak season
    ('Chañaral coast', 'Chile', -26.6, -70.71, 290, (-26.5, -71.0)),  # 58% year-round, 65% peak season
    ('Greymouth coast', 'New Zealand', -42.5, 171.19, 270, (-42.5, 171.0)),  # 52% year-round, 59% peak season
    ('Marcona coast', 'Peru', -15.7, -74.55, 208, (-16.5, -75.0)),  # 51% year-round, 65% peak season
    ('Port Lincoln coast', 'Australia', -34.54, 135.13, 289, (-34.5, 135.0)),  # 47% year-round, 54% peak season
    ('Meningie coast', 'Australia', -36.1, 139.55, 227, (-36.5, 139.0)),  # 47% year-round, 57% peak season
    ('Oranjemund coast', 'Namibia', -28.54, 16.34, 277, (-28.5, 16.0)),  # 46% year-round, 55% peak season
    ('Nzeto coast', 'Angola', -7.23, 12.86, 253, (-7.5, 12.0)),  # 45% year-round, 67% peak season
    ('Quellón coast', 'Chile', -43.57, -74.84, 274, (-43.5, -76.0)),  # 45% year-round, 58% peak season
    ('Port Denison coast', 'Australia', -29.14, 114.86, 244, (-29.5, 114.0)),  # 43% year-round, 50% peak season
    ('Blitar coast', 'Indonesia', -8.24, 111.06, 192, (-8.5, 111.0)),  # 42% year-round, 58% peak season
    ('Nelson coast', 'New Zealand', -40.95, 172.14, 297, (-40.5, 171.0)),  # 41% year-round, 50% peak season
    ('Kaitaia coast', 'New Zealand', -34.46, 172.71, 265, (-34.5, 172.0)),  # 41% year-round, 48% peak season
    ('Bredasdorp coast', 'South Africa', -34.44, 20.77, 170, (-35.5, 21.0)),  # 40% year-round, 46% peak season
    ('Pointe-Noire coast', 'Congo (Brazzaville)', -4.2, 11.39, 233, (-4.5, 11.0)),  # 40% year-round, 62% peak season
    ('Huarmey coast', 'Peru', -10.09, -78.19, 243, (-10.5, -79.0)),  # 40% year-round, 52% peak season
    ('Gwadar coast', 'Pakistan', 25.25, 63.17, 192, (24.5, 63.0)),  # 24% year-round, 82% peak season
    ('Salalah coast', 'Oman', 16.72, 53.3, 232, (16.5, 53.0)),  # 21% year-round, 72% peak season
    ('Karachi coast', 'Pakistan', 24.93, 66.68, 235, (24.5, 66.0)),  # 19% year-round, 66% peak season
    ('Chabahar coast', 'Iran', 25.14, 61.24, 199, (24.5, 61.0)),  # 19% year-round, 69% peak season
]
