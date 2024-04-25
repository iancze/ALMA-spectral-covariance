import numpy as np
import matplotlib.pyplot as plt 

# demonstrate that we can properly Hann-smooth a sequence in python.

N = 16

vals = np.zeros(N)

vals[N//2] = 1.

hann_kernel = np.array([0.25, 0.5, 0.25])

result = np.convolve(vals, hann_kernel, mode="same")


fig, ax = plt.subplots(nrows=2, sharex=True)
ax[0].plot(vals, "o")
ax[1].plot(result, "o")
fig.savefig("convolution-test.png", dpi=200)
