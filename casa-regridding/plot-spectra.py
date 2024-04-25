import mpoldatasets.image
import matplotlib.pyplot as plt
import numpy as np

# load the FITS file
# extract the peak pixel value for each channel
# plot

def get_spectrum(fname):

    data, RA, DEC, ext, vels, beam = mpoldatasets.image.load_cube(fname)

    # flip the data back
    data = data[::-1]
    vels = vels[::-1]

    spectrum = np.max(data, axis=(1,2))
    return spectrum 


fig, ax = plt.subplots(nrows=1)

ax.plot(get_spectrum("topo-wavelet-linear-t0.fits"), "o", label="linear topo")
ax.plot(get_spectrum("bary-wavelet-linear.fits"), "o", label="linear")
ax.plot(get_spectrum("bary-wavelet-cubic.fits"), "o", label="cubic")
ax.legend()
fig.savefig("spectrum-bary.png", dpi=200)
