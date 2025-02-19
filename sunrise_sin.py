from numpy import sin, abs, array, max, pi, sin, cos, acos, append

from scipy.optimize import curve_fit

import json

with open('comparison.json', 'r') as fp:
    comparison = json.load(fp)

def deg2rad(d):
    return d / 180 * pi

class Parameters:
    def __init__(self, latitude_deg, longitude_deg):
        self.longitude_deg = longitude_deg
        self.longitude_days = self.longitude_deg / 360

parameters = Parameters(60.1708, 24.9375)

# i = 0 ~= 1950-01-01T12:00
def i2jd(i):
    return i + 2433283

def jd2i(jd):
    return jd - 2433283

def i2readable(i):
    if i < 0:
        return '-' + u2readable(-i)
    else:
        return u2readable(i)

def u2readable(i):
    d = int(i)
    r = i - d
    s = ''
    if d:
        s += '%d d ' % d
    r *= 24
    h = int(r)
    r -= h
    if h:
        s += '%d h ' % h
    r *= 60
    m = int(r)
    r -= m
    if m:
        s += '%d m ' % m
    r *= 60
    if r:
        s += '%s s' % r
    return s

def sin_daylight(x, a, b, c, d):
    return c * sin(a * (x - b)) + d

def sunchange(jd, is_sunset):
    i = jd2i(jd)
    daylight = sin_daylight(i, *sin_params)
    midday = int(jd + 0.5) - parameters.longitude_days
    if is_sunset:
        return midday + daylight / 2
    else:
        return midday - daylight / 2

def sunrise(jd):
    return sunchange(jd, False)

def sunset(jd):
    return sunchange(jd, True)

def measure_error():
    errors = []
    for i, (proto_sunrise, proto_sunset) in enumerate(comparison):
        if i > 27375:
            break
        jd = i2jd(i)
        our_sunrise = sunrise(jd)
        our_sunset = sunset(jd)
        error = our_sunrise - proto_sunrise
        errors.append(error)
    print('min error %s' % i2readable(min(errors)))
    print('max error %s' % i2readable(max(errors)))

def fit_params():
    domain = array(range(len(comparison)))
    v = array([item[1] - item[0] for item in comparison])
    init_params = [0.0172142, 80, 0.25, 0.5]
    sin_params = curve_fit(sin_daylight, domain, v, init_params)[0]
    print(sin_params)
    return sin_params

sin_params = fit_params()
measure_error()

