import casatools
import numpy as np
import matplotlib.pyplot as plt

# initialize tools
tb = casatools.table()
ms = casatools.ms()
msmd = casatools.msmetadata()

fname = "bary-wavelet.ms"

spw_id = 0
msmd.open(fname)
chan_freq = msmd.chanfreqs(spw_id)
# get the direction for a field
msmd.done()

ms.open(fname, nomodify=False)
ms.selectinit(datadescid=spw_id)
ms.select({"scan_number": [106, 106]})
query = ms.getdata(["data"], ifraxis=True)

data = query["data"]
print(data.shape)
# 2, 3840, 48165
# npol, nchan, nvis

ms.close()

# (2, 3840, 741, 65)
npol, nchan, nbaselines, ntimes = data.shape

fig, ax = plt.subplots(nrows=3, sharex=True)
for i in range(5):
    ax[0].plot(chan_freq, data[0,:,0,i].real, ".")
    ax[1].plot(chan_freq, data[0,:,0,i].imag, ".")
    ax[2].plot(chan_freq, np.abs(data[0,:,0,i]), ".")
ax[2].set_xlim(345.7535e9, 345.754e9)
ax[0].set_ylabel(r"$\Re f$")
ax[1].set_ylabel(r"$\Im f$")
ax[2].set_ylabel(r"$|f|$")
ax[2].set_xlabel(r"$t$")
fig.savefig("wavelet-function-bary-t0.png", dpi=200)
