import casatasks
import mpoldatasets.image

# Just create a simple image, and then plot 
# should span about 10 channels
# print("central channel", t0/Dt + N/2) # 2021

# LkCa 15 imaging parameters 
# restfreq = '345.79598990GHz'
# extent = 10 # km/s
# vsys = 6.3 # km/s                    # LSR systemic velocity in km/s

# start = '{:.3f} km/s'.format(-extent + vsys)
# width = 0.028 # km/s
# nchan = int((2.0 * extent) / width)

# nu0 of the line = nu0 345765642321.9087
nu0 = 345765642321.9087
start_nu = nu0 - (15 * 15.259)

vis = "bary-wavelet.ms"
imagename = "bary-wavelet-linear"

mpoldatasets.image.clear_extensions(imagename)

# image parameters

casatasks.tclean(
    vis=vis,
    imagename=imagename,
    specmode="cube",
    weighting="briggs",
    robust=0.0,
    imsize=1024,
    cell=".02arcsec",
    niter=0,
    threshold="0.05mJy",
    nchan=50,
    start=2780,
    width=1,
    veltype="radio",
    outframe="bary",
    interpolation="linear"
    # timerange='17:19:46~17:20:30' # single integration
)

mpoldatasets.image.exportfits(imagename + ".image", imagename + ".fits")
