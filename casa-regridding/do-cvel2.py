import casatasks
import casatools
import os
import shutil

# perform cvel2 of bary-wavelet.ms to BARY frame, for linear, and cubic interpolation.

fname = "bary-wavelet.ms"

tb = casatools.table()
ms = casatools.ms()

msmd = casatools.msmetadata()

for window in ["rect", "hann"]:

    # get a list of spectral windows 
    for interp in ["linear", "cubic"]:
        outvis = "{:}-{:}.ms".format(window, interp)

        # remove outvis, if exists
        if os.path.exists(outvis):
            shutil.rmtree(outvis)

        fname = "bary-wavelet-{:}.ms".format(window)
        casatasks.cvel2(fname, outputvis=outvis, outframe="BARY", interpolation=interp, mode="frequency")
