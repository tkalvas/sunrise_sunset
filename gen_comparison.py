from skyfield import almanac
from skyfield.api import N, S, E, W, load, wgs84

import json

ts = load.timescale()
eph = load('de421.bsp')
sun = eph['Sun']
helsinki = wgs84.latlon(60.1708 * N, 24.9375 * E)
observer = eph['Earth'] + helsinki

t0 = ts.utc(1950, 1, 1, -2)
t1 = ts.utc(2050, 1, 1, -2)
t, y = almanac.find_risings(observer, sun, t0, t1)
u, z = almanac.find_settings(observer, sun, t0, t1)

if len(t) != len(u):
    raise Exception('|t| != |u|')

rows = [[t[i].tdb, u[i].tdb] for i in range(len(t))]

with open('comparison.json', 'w') as fp:
    json.dump(rows, fp)
