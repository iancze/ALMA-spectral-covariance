import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sici

def sinc_SRF(x):
    return np.sinc(x)

def sinc_averaged_SRF(x, navg=2):
    si1, _ = sici(np.pi * (x + navg/2))
    si2, _ = sici(np.pi * (x - navg/2))

    return (si1 - si2) / np.pi

def hann_SRF(x):
    return 2 * (0.5 * np.sinc(x) + 0.25 * np.sinc(x + 1) + 0.25 * np.sinc(x - 1))

def hann_averaged_SRF(x, navg=2):
    si1, _ = sici(np.pi * (x + navg/2))
    si2, _ = sici(np.pi * (x - navg/2))
    term1 = 0.5/np.pi * (si1 - si2)

    si3, _ = sici(np.pi * (x + 1 + navg/2))
    si4, _ = sici(np.pi * (x + 1 - navg/2))
    term2 = 0.25/np.pi * (si3 - si4)

    si5, _ = sici(np.pi * (x - 1 + navg/2))
    si6, _ = sici(np.pi * (x - 1 - navg/2))
    term3 = 0.25/np.pi * (si5 - si6)

    return term1 + term2 + term3


xs = np.linspace(0, 50, num=1000)
fig, ax = plt.subplots(nrows=2, sharex=True, figsize=(5,4))

ax[0].plot(xs, sinc_SRF(xs), label="rect")
ax[0].plot(xs, hann_SRF(xs), label="Hann")
ax[0].legend()
ax[0].set_xlabel("channel")

for navg in [2, 4, 8]:
    ax[1].plot(xs/navg, sinc_averaged_SRF(xs, navg=navg), color="C0")
    ax[1].plot(xs/navg, hann_averaged_SRF(xs, navg=navg), color="C1")

ax[1].set_xlabel("averaged channel")

ax[1].set_xlim(0, 4)
fig.savefig("SRF.png", dpi=300)
