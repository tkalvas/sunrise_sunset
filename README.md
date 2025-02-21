# sunrise_sunset
Rough tests on sunrise-sunset calculation methods

## Raison d'être

Every now and then someone comes up with the idea that sunrise time follows a sine curve. Well, it depends on what kind of error you are willing to tolerate. But it's not a good approximation, and geometrically baseless.

## Files

### gen_comparison.py

Generates a dump of sunrise and sunset times for each day from 1950 to 2050, using the Skyfield astronomical library for Python.

### sunrise_sin.py

Finds a good parametrization of the sine function for the purpose of sunrise and sunset forecasting using curve_fit from the SciPy library. Calculates the maximum error over the comparison data. (I know, don't use the same data to build a method and test it.)

### sunrise_short.py

A method of calculating the sunrise and sunset times based on a short literature review. Fast, unlike Skyfield.

### sunrise_parametrizable.py

A version of the method in sunrise_short.py, where the individual parts of the algorithm can be toggled on and off separately, for the purpose of finding out the effect on the result of each of them.

## Accuracy of the results

sunrise_sin error ranges from -15:22 to 26:32 (minutes:seconds). sunrise_short error ranges from -7.8 to 5.8 seconds.

The effect of the separate parts of the algorithm on the error varies. Not doing the third sunrise angle round, using a less accurate true anomaly, and ignoring orbital elements changing in time are under a minute. The less accurate equation of time, only one sunrise angle round, and doing no delta T corrections are under five minutes. Using a linear declination calculation, ignoring orbital eccentricity, atmospheric diffraction at the horizon, and the fact that the sun is not a point light source are under ten minutes. Not performing any kind of equation of time calculation gives an error of slightly more than 30 minutes. The details are in sunrise_errors.txt

There is a physical limitation on defining the exact time of sunrise and sunset. It is defined by the upper limb of the sun being visible, which means we have to take into account atmospheric refraction. The actual amount of refraction depends on weather conditions, so seconds in sunrise and sunset times are pretty much meaningless.

The analysis is done on the locations of Helsinki, latitude 60.1708 degrees North, and Punta Arenas, latitude 53.1667 degrees South. There is reason to believe errors of any realistic method get worse the closer to the poles you get. Note that I didn't bother handling the case where sun does not set (or rise) at all, which occasionally happens if you are north of about 66.6 N or south of 66.6 S.

## A note on time scales

All calculations are in TBD, Barycentric Dynamical Time, which is historical and superseded but pretty much the same as TT, Terrestiol Time. TT is not TAI, and TAI is not UTC. If you want to use this code to actually publish sunrise and sunset times, you probably need to convert to UTC, which means dealing with leap seconds, so I left it out.
