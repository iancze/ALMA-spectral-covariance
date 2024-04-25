import casatools
import numpy as np
import wavelets
import barytimes
import visread.utils


# initialize tools
ms = casatools.ms()
msmd = casatools.msmetadata()

fname = "bary-wavelet-rect.ms"

spw_id = 0
msmd.open(fname)
chan_freq = msmd.chanfreqs(spw_id)
msmd.done()

ms.open(fname, nomodify=False)
ms.selectinit(datadescid=spw_id)
query = ms.getdata(["data"], ifraxis=True)

data = query["data"]
print(data.shape)
# 2, 3840, 48165
# npol, nchan, nvis

# (2, 3840, 741, 65)
npol, nchan, nbaselines, ntimes = data.shape

# XX, YY polarizations are a weighted-averaged
# means scale factor should still be the same

# we are putting in a point source, which means the 
# visibility will be the same for all baselines

# BARY frame coordinates of line
N = 3840
x0 = 345.765e9 # Hz
f0 = 1.7e-05 # oscillations per Hz
sigma = 30000 # Hz

# normalize the prefactor. Max of 1 Jy
prefactor = sigma * np.sqrt(2 * np.pi)

# create the array to hold data
data = np.empty(data.shape, dtype=np.complex128)
broadcast = np.ones((2, nchan, nbaselines), dtype=np.complex128)

for i in range(len(barytimes.corrs)):
    # for each timestamp, apply the barycentric correction
    # to go from TOPO to BARY frame, for this epoch
    # barrycorpy convention is to add
    chan_freq_epoch = visread.utils.doppler_shift(chan_freq, barytimes.corrs[i])

    # evaluate wavelet on shifted freqs
    # 1 Jy max
    hs = prefactor * wavelets.raisedcosgauss(chan_freq_epoch, x0, f0, sigma)

    # insert into dataset along time axis
    data[:,:,:,i] = hs[np.newaxis, :, np.newaxis] * broadcast

# collapse the ifr axis
# refold data
data = np.reshape(data, (2, nchan, nbaselines * ntimes), order="F")

ms.putdata({"data":data})
ms.close()
