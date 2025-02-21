import json
from math import sin, cos, tan, asin, acos, atan, pi, atan2

with open('comparison.json', 'r') as fp:
    comparison = json.load(fp)

def deg2rad(d):
    return d / 180 * pi

def crop(rad):
    return rad % (2*pi)

def crop2(rad):
    return crop(rad + pi) - pi

def jd2t(jd):
    return (jd - 2451545) / 365250

NONE = 0
SOME = 1
PROPER = 2

class AlgorithmOptions:
    def __init__(self, sunrise_angle_rounds=3, sun_radius=True,
                 horizon_diffraction=True, delta_t=True,
                 orbital_elements=True, orbital_eccentricity=True,
                 equation_of_time=PROPER, nonlinear_declination=True,
                 accurate_true_anomaly=True):
        self.sunrise_angle_rounds = sunrise_angle_rounds
        self.sun_radius = sun_radius
        self.horizon_diffraction = horizon_diffraction
        self.delta_t = delta_t
        self.orbital_elements = orbital_elements
        self.orbital_eccentricity = orbital_eccentricity
        self.equation_of_time = equation_of_time
        self.nonlinear_declination = nonlinear_declination
        self.accurate_true_anomaly = accurate_true_anomaly

class Parameters:
    def __init__(self, latitude_deg, longitude_deg):
        self.latitude_deg = latitude_deg
        self.longitude_deg = longitude_deg
        self.latitude_rad_prime = pi/2 - deg2rad(self.latitude_deg)
        self.longitude_days = self.longitude_deg / 360

class DegElement:
    def __init__(self, deg, *secs):
        self.deg = deg
        self.secs = secs

    def deg_at(self, t):
        tt = 1.0
        acc = 0.0
        for sec in self.secs:
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
        tt = 1.0
        acc = 0.0
        for corr in self.corrs:
            tt *= t
            acc += tt * corr
        return self.base + acc

class EarthOrbit:
    # from https://articles.adsabs.harvard.edu/pdf/1994A%26A...282..663S
    # t is TDB in thousands of julian years from J2000
    mean_longitude = DegElement(100.46645683, 1296027711.03429, 109.15809,
                                0.07207, -0.23530, -0.00180, 0.00020)
    eccentricity = Element(0.0167086342, -0.0004203654, -0.0000126734,
                           1444e-10, -2e-10, 3e-10)
    longitude_of_perihelion = DegElement(102.93734808, 61900.55290,
                                         164.47797, -0.06365, -0.12090,
                                         0.00298, 0.00020)
    obliquity_of_the_ecliptic_rad = 0.4090926294954782

class FixedEarthOrbit:
    mean_longitude = DegElement(100.46645683, 1296027711.03429)
    eccentricity = Element(0.0167086342)
    longitude_of_perihelion = DegElement(102.93734808)
    obliquity_of_the_ecliptic_rad = 0.4090926294954782

earth_orbit = EarthOrbit()
fixed_earth_orbit = FixedEarthOrbit()

class OrbitPosition:
    def __init__(self, algorithm, jd):
        self.jd = jd
        self.t = jd2t(jd)
        self.mean_longitude = algorithm.mean_longitude_rad(self.t)
        self.eccentricity = algorithm.eccentricity(self.t)
        self.longitude_of_perihelion = algorithm.longitude_of_perihelion_rad(self.t)
        self.mean_anomaly = crop(self.mean_longitude - self.longitude_of_perihelion)
        self.true_anomaly = algorithm.true_anomaly(self.mean_anomaly, self.eccentricity)
        self.true_longitude = crop(self.longitude_of_perihelion + self.true_anomaly)

class Algorithm:
    def __init__(self, options, parameters):
        self.options = options
        self.parameters = parameters

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

    def delta_t_seconds(self, t):
        if self.options.delta_t:
            y = 1000.0 * t + 2000
            delta_t_recent = self.delta_t_recent
            for i in range(len(delta_t_recent) - 1):
                if delta_t_recent[i][0] <= y < delta_t_recent[i + 1][0]:
                    x = y - delta_t_recent[i][0]
                    k = (delta_t_recent[i + 1][1] - delta_t_recent[i][1]) / 5.0
                    return delta_t_recent[i][1] + k * x
            raise Exception('input not supported ' + str(t))
        else:
            return 0

    def delta_t_days(self, t):
        return self.delta_t_seconds(t) / 86400.0

    def sunrise_depth_rad(self):
        val = pi/2
        if self.options.sun_radius:
            val += deg2rad(34.0/60.0)
        if self.options.horizon_diffraction:
            val += deg2rad(16.0/60.0)
        return val

    def obliquity_of_the_ecliptic_rad(self, t):
        if self.options.orbital_elements:
            return earth_orbit.obliquity_of_the_ecliptic_rad
        else:
            return fixed_earth_orbit.obliquity_of_the_ecliptic_rad

    def sun_declination(self, t, ecliptic_longitude):
        if self.options.nonlinear_declination:
            return asin(sin(self.obliquity_of_the_ecliptic_rad(t)) * sin(ecliptic_longitude))
        else:
            return self.obliquity_of_the_ecliptic_rad(t) * sin(ecliptic_longitude)

    def sunrise_angle(self, t, ecliptic_longitude):
        latP = self.parameters.latitude_rad_prime
        declination = self.sun_declination(t, ecliptic_longitude)
        depth = self.sunrise_depth_rad()
        return acos((cos(depth) + cos(latP) * sin(declination)) / (sin(latP) * cos(declination)))

    def mean_longitude_rad(self, t):
        if self.options.orbital_elements:
            return earth_orbit.mean_longitude.rad_at(t)
        else:
            return fixed_earth_orbit.mean_longitude.rad_at(t)

    def eccentricity(self, t):
        if not self.options.orbital_eccentricity:
            return 0.0
        if self.options.orbital_elements:
            return earth_orbit.eccentricity.at(t)
        else:
            return fixed_earth_orbit.eccentricity.at(t)

    def longitude_of_perihelion_rad(self, t):
        if self.options.orbital_elements:
            return earth_orbit.longitude_of_perihelion.rad_at(t)
        else:
            return fixed_earth_orbit.longitude_of_perihelion.rad_at(t)

    def right_ascension(self, t, lon):
        # avoid quadrant mess with atan(cos(eps)*tan(lon))
        x = sin(lon)
        y = cos(lon)
        return atan2(cos(self.obliquity_of_the_ecliptic_rad(t)) * x, y)

    def true_anomaly(self, m, e):
        if self.options.accurate_true_anomaly:
            return (m + (2 * e - e**3 / 4) * sin(m) +
                    5 * e**2 / 4 * sin(2 * m) +
                    13 * e**3 / 12 * sin(3 * m))
        else:
            return m + 2 * e * sin(m)

    def equation_of_time(self, position):
        if self.options.equation_of_time == PROPER:
            return crop2(position.mean_anomaly + position.longitude_of_perihelion - self.right_ascension(position.t, position.true_longitude))
        elif self.options.equation_of_time == SOME:
            d = 6.24004077 + 0.01720197 * (position.jd - 2451545)
            mins = -7.659 * sin(d) + 9.863 * sin(2*d + 3.5932)
            return mins / 60 / 24 * (2*pi)
        else:
            return 0.0

    def equation_of_time_days(self, position):
        return self.equation_of_time(position) / (2*pi)

    def sunchange(self, jd, is_set):
        raw_midday = int(jd + 0.5)
        for i in range(self.options.sunrise_angle_rounds):
            position = OrbitPosition(self, jd)
            dt = self.delta_t_days(position.t)
            eot = self.equation_of_time_days(position)
            half_daylength = self.sunrise_angle(position.t, position.true_longitude) / (2*pi)
            jd = raw_midday - self.parameters.longitude_days - eot + dt
            if is_set:
                jd += half_daylength
            else:
                jd -= half_daylength
        return jd

    def sunrise(self, jd):
        return self.sunchange(jd, False)

    def sunset(self, jd):
        return self.sunchange(jd, True)

### only used for the error measurement:

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

# i = 0 ~= 1950-01-01T12:00
def i2jd(i):
    return i + 2433283

def measure_error():
    parameters = Parameters(60.1708, 24.9375)
    #parameters = Parameters(-53.1667, -70.9333)
    algorithm_options = AlgorithmOptions()
    algorithm = Algorithm(algorithm_options, parameters)
    errors = []
    for i, (proto_sunrise, proto_sunset) in enumerate(comparison):
        # we don't know delta T of the future yet
        if i > 27375:
            break
        jd = i2jd(i)
        our_sunrise = algorithm.sunrise(jd)
        our_sunset = algorithm.sunset(jd)
        error = our_sunrise - proto_sunrise
        errors.append(error)
    print('min error %s' % i2readable(min(errors)))
    print('max error %s' % i2readable(max(errors)))

if __name__ == '__main__':
    measure_error()
