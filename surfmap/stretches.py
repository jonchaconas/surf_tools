"""Under-the-radar surf coast found by discover.py (regenerate with `python -m surfmap discover`).

Each entry is (name, region, lat, lon, facing, cell): a coast point, the compass direction
it faces out to sea, and the WW3 1 degree cell its figures come from. Numbering on the
maps follows this order (best year-round first).
"""

STRETCHES = [
    ('Ancud coast', 'Chile', -41.45, -73.86, 244, (-41.5, -74.0)),  # 80% year-round, 87% peak season
    ('Camaná coast', 'Peru', -16.61, -72.79, 193, (-17.5, -73.0)),  # 79% year-round, 85% peak season
    ('Antofagasta coast', 'Chile', -24.33, -70.55, 248, (-24.5, -71.0)),  # 79% year-round, 84% peak season
    ('La Ligua coast', 'Chile', -31.5, -71.58, 269, (-31.5, -72.0)),  # 73% year-round, 78% peak season
    ('Toltén coast', 'Chile', -38.51, -73.52, 271, (-38.5, -74.0)),  # 70% year-round, 75% peak season
    ('Isabela coast', 'Galápagos, Ecuador', -0.39, -91.65, 253, (-0.5, -92.0)),  # 69% year-round, 80% peak season
    ('Chañaral coast', 'Chile', -26.6, -70.71, 290, (-26.5, -71.0)),  # 67% year-round, 71% peak season
    ('Soyo coast', 'Angola', -6.19, 12.33, 227, (-6.5, 12.0)),  # 64% year-round, 77% peak season
    ('Sumbe coast', 'Angola', -10.26, 13.5, 243, (-10.5, 13.0)),  # 64% year-round, 77% peak season
    ('Pointe-Noire coast', 'Congo (Brazzaville)', -4.2, 11.39, 233, (-4.5, 11.0)),  # 63% year-round, 78% peak season
    ('Greymouth coast', 'New Zealand', -42.5, 171.19, 270, (-42.5, 171.0)),  # 60% year-round, 66% peak season
    ('La Serena coast', 'Chile', -29.35, -71.35, 255, (-29.5, -72.0)),  # 59% year-round, 66% peak season
    ('Omboué coast', 'Gabon', -2.03, 9.4, 220, (-2.5, 9.0)),  # 57% year-round, 73% peak season
    ('Escuintla coast', 'Mexico', 15.14, -92.81, 196, (14.5, -93.0)),  # 52% year-round, 67% peak season
    ('Marcona coast', 'Peru', -15.7, -74.55, 208, (-16.5, -75.0)),  # 51% year-round, 69% peak season
    ('Acapulco coast', 'Mexico', 16.92, -100.02, 177, (16.5, -100.0)),  # 51% year-round, 64% peak season
    ('Port Lincoln coast', 'Australia', -34.54, 135.13, 289, (-34.5, 135.0)),  # 50% year-round, 54% peak season
    ('Luanda coast', 'Angola', -8.55, 13.37, 279, (-8.5, 13.0)),  # 50% year-round, 73% peak season
    ('Meningie coast', 'Australia', -36.1, 139.55, 227, (-36.5, 139.0)),  # 49% year-round, 59% peak season
    ('Oranjemund coast', 'Namibia', -28.54, 16.34, 277, (-28.5, 16.0)),  # 48% year-round, 58% peak season
    ('Blitar coast', 'Indonesia', -8.24, 111.06, 192, (-8.5, 111.0)),  # 48% year-round, 65% peak season
    ('Maintirano coast', 'Madagascar', -18.5, 44.11, 272, (-18.5, 44.0)),  # 46% year-round, 55% peak season
    ('Morombe coast', 'Madagascar', -22.38, 43.26, 245, (-22.5, 43.0)),  # 45% year-round, 55% peak season
    ('São Tomé coast', 'Sao Tome and Principe', 0.07, 6.52, 223, (-0.5, 6.0)),  # 45% year-round, 58% peak season
    ('Huarmey coast', 'Peru', -10.09, -78.19, 243, (-10.5, -79.0)),  # 44% year-round, 58% peak season
    ('Quellón coast', 'Chile', -43.57, -74.84, 274, (-43.5, -76.0)),  # 43% year-round, 56% peak season
    ('Rabat coast', 'Morocco', 34.13, -6.76, 331, (34.5, -7.0)),  # 43% year-round, 65% peak season
    ('Juradó coast', 'Colombia', 6.69, -77.53, 248, (6.5, -78.0)),  # 42% year-round, 57% peak season
    ('Port Denison coast', 'Australia', -29.14, 114.86, 244, (-29.5, 114.0)),  # 42% year-round, 50% peak season
    ('Kaitaia coast', 'New Zealand', -34.46, 172.71, 265, (-34.5, 172.0)),  # 42% year-round, 49% peak season
    ('Nelson coast', 'New Zealand', -40.95, 172.14, 297, (-40.5, 171.0)),  # 42% year-round, 52% peak season
    ('San Cristóbal coast', 'Galápagos, Ecuador', -0.91, -89.42, 145, (-1.5, -89.0)),  # 41% year-round, 66% peak season
    ('Kangaroo Island coast', 'Australia', -35.94, 136.59, 220, (-36.5, 136.0)),  # 40% year-round, 48% peak season
]
