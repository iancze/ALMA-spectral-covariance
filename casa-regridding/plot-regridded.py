import casatools
import numpy as np
import wavelets
import matplotlib.pyplot as plt

# initialize tools
tb = casatools.table()
ms = casatools.ms()
msmd = casatools.msmetadata()

# BARY frame coordinates of line
N = 3840
x0 = 345.765e9 # Hz
f0 = 1.7e-05 # oscillations per Hz
sigma = 30000 # Hz

# normalize the prefactor. Max of 2 Jy
prefactor = sigma * np.sqrt(2 * np.pi)


for window in ["rect", "hann"]:

    for interp in ["linear", "cubic"]:

        fname = "{:}-{:}.ms".format(window, interp)
        spw_id = 0
        msmd.open(fname)
        chan_freq = msmd.chanfreqs(spw_id)
        # get the direction for a field
        msmd.done()

        ms.open(fname)
        ms.selectinit(datadescid=spw_id)
        # ms.select({"scan_number": [106, 106]})
        query = ms.getdata(["data"], ifraxis=True)

        data = query["data"]
        print(data.shape)

        ms.close()

        # (2, 3840, 741, 65)
        npol, nchan, nbaselines, ntimes = data.shape

        xs_fine = np.linspace(np.min(chan_freq), np.max(chan_freq), num=100000)
        hs_fine = prefactor * wavelets.raisedcosgauss(xs_fine, x0, f0, sigma)

        hs = prefactor * wavelets.raisedcosgauss(chan_freq, x0, f0, sigma)
        hs_hann = wavelets.apodize(hs)
        xs_hann_fine, hs_hann_fine = wavelets.apodize_fine(chan_freq, hs)

        fig, ax = plt.subplots(nrows=3, sharex=True)

        if window == "rect":
            ax[0].plot(xs_fine, hs_fine, label="rect true")
        if window == "hann":
            ax[0].plot(xs_hann_fine, hs_hann_fine.real, label="hann true")
        
        for i in range(ntimes):
            ax[0].plot(chan_freq, data[0,:,0,i].real, ".", color="k", label="data")
            ax[1].plot(chan_freq, data[0,:,0,i].imag, ".", color="k", label="data")

        # get the average values, too
        data_avg = np.average(data, axis=(0,2,3))
        ax[0].plot(chan_freq, data_avg.real, "o", label="avg")
        ax[1].plot(chan_freq, data_avg.imag, "o", label="avg")

        if window == "rect":
            ax[2].plot(chan_freq, (hs - data_avg).real, "o", label="avg")
        if window == "hann":
            ax[2].plot(chan_freq, (hs_hann - data_avg).real, "o", label="avg")

        ax[2].set_xlim(x0 - 10 * sigma, x0 + 10 * sigma)
        ax[0].set_ylabel(r"$\Re f$")
        ax[1].set_ylabel(r"$\Im f$")
        ax[2].set_ylabel(r"$|f|$")
        ax[2].set_xlabel(r"$\nu$ [BARY]")
        fig.savefig("interpolated-{:}-{:}.png".format(window, interp), dpi=200)
