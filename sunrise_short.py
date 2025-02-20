import json
from math import sin, cos, tan, asin, acos, atan, pi, atan2
from types import SimpleNamespace

with open('comparison.json', 'r') as fp:
    comparison = json.load(fp)

def rad2deg(r):
    return 180 * r / pi

def deg2rad(d):
    return d / 180 * pi

class Parameters:
    def __init__(self, latitude_deg, longitude_deg):
        self.latitude_deg = latitude_deg
        self.longitude_deg = longitude_deg
        self.latitude_rad_prime = pi/2 - deg2rad(self.latitude_deg)
        self.longitude_days = self.longitude_deg / 360
        self.obliquity_of_the_ecliptic_rad = 0.4090926294954782
        self.sunrise_depth_rad = 1.5853407372281827 # 90 deg 50 min

parameters = Parameters(60.1708, 24.9375)

def jd2t(jd):
    return (jd - 2451545) / 365250

# i = 0 ~= 1950-01-01T12:00
def i2jd(i):
    return i + 2433283

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

delta_t_recent = [
    [1950, 29],
    [1955, 31.1],
    [1960, 33.2],
    [1965, 35.7],
    [1970, 40.2],
    [1975, 45.5],
    [1980, 50.5],
    [1985, 54.3],
    [1990, 56.9],
    [1995, 60.8],
    [2000, 63.8],
    [2005, 64.7],
    [2010, 66.1],
    [2015, 67.7],
    [2020, 69.3],
    [2025, 69.1]
]

def delta_t_seconds(t):
    y = 1000.0 * t + 2000
    for i in range(len(delta_t_recent) - 1):
        if delta_t_recent[i][0] <= y < delta_t_recent[i + 1][0]:
            x = y - delta_t_recent[i][0]
            k = (delta_t_recent[i + 1][1] - delta_t_recent[i][1]) / 5.0
            return delta_t_recent[i][1] + k * x
    raise Exception('input not supported ' + str(t))

def delta_t_days(t):
    return delta_t_seconds(t) / 86400.0

class DegElement:
    def __init__(self, deg, *secs):
        self.deg = deg
        self.secs = secs

    def deg_at(self, t):
        tt = t
        acc = self.secs[0] * t
        for sec in self.secs[1:]:
            tt *= t
            acc += tt * sec
        return self.deg + acc / 3600

    def rad_at(self, t):
        return deg2rad(self.deg_at(t))

class Element:
    def __init__(self, base, *corrs):
        self.base = base
        self.corrs = corrs

    def at(self, t):
        tt = t
        acc = self.corrs[0] * t
        for corr in self.corrs[1:]:
            tt *= t
            acc += tt * corr
        return self.base + acc

class EarthOrbit:
    # from https://articles.adsabs.harvard.edu/pdf/1994A%26A...282..663S
    # t is TDB in thousands of julian years from J2000
    semi_major_axis = 1.00000178
    mean_longitude = DegElement(100.46645683, 1296027711.03429, 109.15809,
                                0.07207, -0.23530, -0.00180, 0.00020)
    eccentricity = Element(0.0167086342, -0.0004203654, -0.0000126734,
                           1444e-10, -2e-10, 3e-10)
    longitude_of_perihelion = DegElement(102.93734808, 61900.55290,
                                         164.47797, -0.06365, -0.12090,
                                         0.00298, 0.00020)

earth_orbit = EarthOrbit()

def crop(rad):
    return rad % (2*pi)

def crop2(rad):
    return crop(rad + pi) - pi

def true_anomaly(m, e):
    return (m + (2 * e - e**3 / 4) * sin(m) +
            5 * e**2 / 4 * sin(2 * m) +
            13 * e**3 / 12 * sin(3 * m))

def props(t):
    ml = earth_orbit.mean_longitude.rad_at(t)
    e = earth_orbit.eccentricity.at(t)
    lp = earth_orbit.longitude_of_perihelion.rad_at(t)
    ma = crop(ml - lp)
    ta = true_anomaly(ma, e)
    tl = crop(lp + ta)
    return SimpleNamespace(t=t, ml=ml, e=e, lp=lp, ma=ma, ta=ta, tl=tl)

def right_ascension(lon):
    # avoid quadrant mess with atan(cos(eps)*tan(lon))
    x = sin(lon)
    y = cos(lon)
    return atan2(cos(parameters.obliquity_of_the_ecliptic_rad) * x, y)

def equation_of_time(lp, ma, tl):
    return crop2(ma + lp - right_ascension(tl))

def equation_of_time_days(lp, ma, tl):
    eot = equation_of_time(lp, ma, tl) / (2*pi)
    #print(ma, lp, tl, eot)
    return eot

# P-(l')--Z    P pole output  Z zenith
#  \     /       l' 90-latitude
#   s   d      s season input 90+axial tilt*sin(year fraction)
#    \ /           d depth 90+diffraction+sun radius
#     S        S sun

# cos d = cos l' cos s + sin l' sin s cos P

# sin l' sin s cos P = cos d - cos l' cos s
# cos P = (cos d - cos l' cos s) / (sin l' sin s)

def spherical_daylight(ecliptic_longitude):
    latP = parameters.latitude_rad_prime
    declination = asin(sin(parameters.obliquity_of_the_ecliptic_rad) * sin(ecliptic_longitude))
    cosP = (cos(parameters.sunrise_depth_rad) + cos(latP) * sin(declination)) / (sin(latP) * cos(declination))
    return acos(cosP) / pi

def sunchange(jd, is_set):
    raw_midday = int(jd + 0.5)
    for i in range(3):
        t = jd2t(jd)
        p = props(t)
        dt = delta_t_days(t)
        eot = equation_of_time_days(p.lp, p.ma, p.tl)
        half_daylength = spherical_daylight(p.tl) / 2
        longitude_days = 24.9375 / 360
        jd = raw_midday - parameters.longitude_days - eot + dt
        if is_set:
            jd += half_daylength
        else:
            jd -= half_daylength
    return jd

def sunrise(jd):
    return sunchange(jd, False)

def sunset(jd):
    return sunchange(jd, True)
    
def measure_error():
    errors = []
    for i, (proto_sunrise, proto_sunset) in enumerate(comparison):
        # we don't know delta T of the future yet
        if i > 27375:
            break
        jd = i2jd(i)
        our_sunrise = sunrise(jd)
        our_sunset = sunset(jd)
        error = our_sunrise - proto_sunrise
        errors.append(error)
    print('min error %s' % i2readable(min(errors)))
    print('max error %s' % i2readable(max(errors)))

if __name__ == '__main__':
    measure_error()
