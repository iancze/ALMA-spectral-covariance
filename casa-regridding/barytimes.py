import casatools
import numpy as np
import astropy.time
import barycorrpy

# initialize tools
ms = casatools.ms()
msmd = casatools.msmetadata()

fname = "LkCa15-spw3.ms"

spw_id = 0
msmd.open(fname)

# get the direction for a field
phasecenter = msmd.phasecenter() # units of radians in J2000
# convert phasecenter to degrees
ra = phasecenter["m0"]["value"] * 180/np.pi # degrees
dec = phasecenter["m1"]["value"] * 180/np.pi # degrees

# get the timestamps for the observations
# these are time centroids
# they are MJD but in units of seconds, rather than days
times = msmd.timesforspws(spw_id)
msmd.done()

# convert from seconds to days
tMJDs = times / (60 * 60 * 24)
ts = astropy.time.Time(tMJDs, format="mjd")

# use barycorrpy to get the velocity corrections for all times 
result = barycorrpy.get_BC_vel(ts, obsname="ALMA", ra=ra, dec=dec)
corrs = result[0] # in m/s 

if __name__ == "__main__":
    print(corrs)
