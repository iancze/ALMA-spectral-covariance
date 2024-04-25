import casatools
import numpy as np
import barytimes
import visread.utils
from astropy.constants import c
import wavelets

# perform cvel2 of bary-wavelet.ms to BARY frame, for linear, and cubic interpolation.

tb = casatools.table()
ms = casatools.ms()

msmd = casatools.msmetadata()

spw_id = 0

fname = "bary-wavelet-rect.ms"
msmd.open(fname)
chan_freq_TOPO = msmd.chanfreqs(spw_id)
msmd.done()

# fname = "rect-linear.ms"
# msmd.open(fname)
# chan_freq_out = msmd.chanfreqs(spw_id)
# msmd.done()

# we have a TOPO frequency grid
# so, shift all bary to that specific TOPO grid

# then, relable TOPO grid with BARY-corr velocities at t0

# and plot relative to true signal

# then, shift all frequencies by an amount to match BARY grid to same CASA grid and swap freq labels

ms.open(fname)
ms.selectinit(datadescid=spw_id)
# ms.select({"scan_number": [106, 106]})
query = ms.getdata(["data"], ifraxis=True)

data = query["data"]
ms.close() 

npol, nchan, nbaselines, ntimes = data.shape
# create the array to hold the shifted data
data_out = np.empty(data.shape, dtype=np.complex128)

chan_freq_BARY0 = visread.utils.doppler_shift(chan_freq_TOPO, barytimes.corrs[0])

for i in range(ntimes):

    chan_freq_BARY = visread.utils.doppler_shift(chan_freq_TOPO, barytimes.corrs[i])
    
    nu_shift = np.average(chan_freq_BARY - chan_freq_BARY0)

    # select the data for this timestamp 
    data_out[:,:,:,i] = wavelets.fftshift_data(chan_freq_BARY, data[:,:,:,i], nu_shift)

np.savez("fftshift-rect.npz", data=data_out, chan_freq=chan_freq_BARY0)