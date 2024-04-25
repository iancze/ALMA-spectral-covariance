import barytimes
import casatools
from astropy.constants import c 

msmd = casatools.msmetadata()

spw_id = 0
fname = "bary-wavelet-rect.ms"
msmd.open(fname)
chan_freq_TOPO = msmd.chanfreqs(spw_id)
msmd.done()

fname = "rect-linear.ms"
msmd.open(fname)
chan_freq_BARY = msmd.chanfreqs(spw_id)
msmd.done()

print(((chan_freq_BARY[0] - chan_freq_TOPO[0]) / chan_freq_BARY[0]) * c.value)
print(barytimes.corrs[0])