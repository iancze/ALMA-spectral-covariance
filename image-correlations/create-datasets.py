import numpy as np 
import matplotlib.pyplot as plt 

d = np.load("snapshot-baselines.npz")
u = d["u"]
v = d["v"]

nbaseline = len(u)

nchan = 16

# we won't be making data with any signal, necessarily. Just noise values.

# assign weight

weight = 2.0 # = 1/sigma^2  [1/Jy^2]
sigma = 1/np.sqrt(weight) # [Jy]

# draw independent noise values for every baseline, for all channels.
V_reals = np.random.normal(scale=sigma, size=(nchan, nbaseline))
V_imags = np.random.normal(scale=sigma, size=(nchan, nbaseline))

# plot an example of a correlated channel
fig, ax = plt.subplots(nrows=2, sharex=True)
ax[0].plot(V_reals[:,0], "o", label="raw")
ax[1].plot(V_imags[:,0], "o", label="raw")

# save dataset
np.savez("uncorrelated_dataset.npz", u=u, v=v, weight=weight, V_reals=V_reals, V_imags=V_imags)

# use np.convolve to convolve across frequency channel with Hann SRF (alternately, draw correlated values)
hann_kernel = np.array([0.25, 0.5, 0.25])

# for each baseline, replace values with correlated values
for i in range(nbaseline):
    V_reals[:, i] = np.convolve(V_reals[:, i], hann_kernel, mode="same")
    V_imags[:, i] = np.convolve(V_imags[:, i], hann_kernel, mode="same")

# modify the weight value to reflect the reduced scatter
weight = 8/3 * weight

# plot an example of a correlated channel
ax[0].plot(V_reals[:,0], ".", label="Hann")
ax[0].legend()
ax[1].plot(V_imags[:,0], ".", label="Hann")
ax[0].set_ylabel(r"$\Re{\mathcal{V}_\nu}$ [Jy]")
ax[1].set_ylabel(r"$\Im{\mathcal{V}_\nu}$ [Jy]")
ax[1].set_xlabel("channel")
fig.savefig("example-channel.png", dpi=300)

# weights need to be modified
# save dataset
np.savez("correlated_dataset.npz", u=u, v=v, weight=weight, V_reals=V_reals, V_imags=V_imags)






