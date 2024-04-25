import numpy as np
import matplotlib.pyplot as plt

def k_rect(v):
    return np.sinc(v)

def k_hann(v):
    return 0.25 * (3/2 * np.sinc(v) + 0.25 * np.sinc(2 - v) + 0.25 * np.sinc(2 + v) + np.sinc(1 - v) + np.sinc(1 + v))

vs = np.linspace(0, 6, num=500)
vs_int = np.arange(7)

fig, ax = plt.subplots(nrows=1)
ax.plot(vs, k_rect(vs), label="rect", color="C0")
ax.plot(vs_int, k_rect(vs_int), "o", color="C0")

ax.plot(vs, k_hann(vs), label="Hann", color="C1")
ax.plot(vs_int, k_hann(vs_int), "o", color="C1")
ax.legend()
ax.set_xlabel("channel")
ax.set_ylabel("covariance")
fig.savefig("ks.png", dpi=300)